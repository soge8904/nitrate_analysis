from pathlib import Path
 
import pandas as pd
import pytest
 
from nitrate_analysis.references import ReferenceData

# tests need to have issues in the loaded excel sheet. Do i make an excel sheet for each issue?

#issues that could arise:
#- no excel file found
#- no summary sheet in excel file
#- can't read excel file --> why?
#

 
# ---------------------------------------------------------------------------
# Helper: write a minimal workbook to a temp path so each test can control
# exactly what "bad" data looks like without touching your real Excel file.
# ---------------------------------------------------------------------------
def make_workbook(tmp_path: Path, rows: list[dict], sheet_name: str = "Reference_summary") -> Path:
    path = tmp_path / "refs.xlsx"
    df = pd.DataFrame(rows)
    df.to_excel(path, sheet_name=sheet_name, index=False)
    return path
 
 
VALID_ROWS = [
    {"name": "mu_H2", "value_eV": -6.77},
    {"name": "mu_H2O", "value_eV": -14.22},
    {"name": "mu_NO", "value_eV": -8.10},
    {"name": "mu_HNO3", "value_eV": -19.5},
    {"name": "dG_rxn_redox", "value_eV": 0.35},
    {"name": "dG_rxn_pka_henry", "value_eV": 0.41},
]
 
 
# ---------------------------------------------------------------------------
# Happy path: valid workbook loads and values come back correctly.
# ---------------------------------------------------------------------------
def test_from_excel_loads_valid_workbook(tmp_path): #pytest calls tmp_path that creates a temp dir when included like here
    #make a workbook with valid rows and load it
    path = make_workbook(tmp_path, VALID_ROWS)
    refs = ReferenceData.from_excel(path)
    #did the function do its job without crashing
 
    assert refs.mu_H2 == pytest.approx(-6.77)
    assert refs.get("mu_H2O") == pytest.approx(-14.22)
    #did the function do its job right
 
# ---------------------------------------------------------------------------
# File-level failures
# ---------------------------------------------------------------------------
def test_missing_file(tmp_path):
    #make different file
    missing_path = tmp_path / "does_not_exist.xlsx"
    #test that it raises right error when loading
    with pytest.raises(FileNotFoundError):
        ReferenceData.from_excel(missing_path)
 
def test_missing_sheet_raises_reference_error(tmp_path):
    path = make_workbook(tmp_path, VALID_ROWS, sheet_name="Some_other_sheet")
    with pytest.raises(ReferenceError):
        ReferenceData.from_excel(path, sheet_name="Reference_summary")
 
# ---------------------------------------------------------------------------
# Column-level failures
# ---------------------------------------------------------------------------
def test_missing_required_column_raises(tmp_path):
    # no "value_eV" column at all
    path = make_workbook(tmp_path, [{"name": "mu_H2"}])
    with pytest.raises(ReferenceError, match="missing required columns"):
        ReferenceData.from_excel(path)
 
# ---------------------------------------------------------------------------
# Row-level failures
# ---------------------------------------------------------------------------
def test_duplicate_names_raise(tmp_path):
    rows = VALID_ROWS + [{"name": "mu_H2", "value_eV": -6.99}]  # duplicate
    path = make_workbook(tmp_path, rows)
    with pytest.raises(ReferenceError, match="Duplicate reference names"):
        ReferenceData.from_excel(path)
 
 
def test_blank_value_raises(tmp_path):
    rows = VALID_ROWS + [{"name": "mu_NO2", "value_eV": None}]
    path = make_workbook(tmp_path, rows)
    with pytest.raises(ReferenceError, match="no numeric value"):
        ReferenceData.from_excel(path)
 
 
def test_non_numeric_value_raises(tmp_path):
    rows = VALID_ROWS + [{"name": "mu_NO2", "value_eV": "not_a_number"}]
    path = make_workbook(tmp_path, rows)
    with pytest.raises(ReferenceError, match="non-numeric value"):
        ReferenceData.from_excel(path)

def test_scientific_notation_value_is_accepted(tmp_path):
    rows = VALID_ROWS + [{"name": "mu_NO2", "value_eV": 1.23e-4}]
    path = make_workbook(tmp_path, rows)
    refs = ReferenceData.from_excel(path)
    assert refs.get("mu_NO2") == pytest.approx(1.23e-4)
 
def test_fully_blank_row_is_silently_dropped(tmp_path):
    # An all-NaN row (e.g. a stray blank Excel row) should NOT raise,
    # because dropna(how="all") removes it before validation runs.
    rows = VALID_ROWS + [{"name": None, "value_eV": None}]
    path = make_workbook(tmp_path, rows)
    refs = ReferenceData.from_excel(path)  # should not raise
    assert len(refs.values) == len(VALID_ROWS)
 
 
# ---------------------------------------------------------------------------
# .get() lookup failures
# ---------------------------------------------------------------------------
def test_get_unknown_name_raises_with_available_list(tmp_path):
    path = make_workbook(tmp_path, VALID_ROWS)
    refs = ReferenceData.from_excel(path)
    with pytest.raises(ReferenceError, match="Available references"):
        refs.get("mu_NOT_A_REAL_SPECIES")
 
 
# ---------------------------------------------------------------------------
# Immutability: does frozen=True actually protect what you think it protects?
# ---------------------------------------------------------------------------
def test_cannot_reassign_values_attribute(tmp_path):
    path = make_workbook(tmp_path, VALID_ROWS)
    refs = ReferenceData.from_excel(path)
    with pytest.raises(Exception):  # dataclasses.FrozenInstanceError
        refs.values = {}
 
 
def test_dict_contents_are_still_mutable_despite_frozen(tmp_path):
    # This is the gotcha flagged earlier: frozen=True stops attribute
    # REASSIGNMENT, but does nothing to stop mutating the dict IN PLACE.
    path = make_workbook(tmp_path, VALID_ROWS)
    refs = ReferenceData.from_excel(path)
    refs.values["mu_H2"] = 9999.0  # this succeeds -- frozen does not block it
    assert refs.mu_H2 == 9999.0
 
 
# ---------------------------------------------------------------------------
# nitrate_cycle_inputs: route selection
# ---------------------------------------------------------------------------
#1 test to check that the wrong route raises an error
# test to assert that the right route returns the right value

def test_nitrate_cycle_inputs_redox_route(tmp_path):
    path = make_workbook(tmp_path, VALID_ROWS)
    refs = ReferenceData.from_excel(path)
    dG = refs.get_nitrate_association_dGrxn("redox")
    assert dG == refs.dG_rxn_redox
 
 
def test_nitrate_cycle_inputs_henry_pka_route(tmp_path):
    path = make_workbook(tmp_path, VALID_ROWS)
    refs = ReferenceData.from_excel(path)
    dG = refs.get_nitrate_association_dGrxn("PKA_henry")  # checks case/whitespace handling
    assert dG == refs.dG_rxn_pka_henry
 
def test_nitrate_cycle_inputs_unknown_route_raises(tmp_path):
    path = make_workbook(tmp_path, VALID_ROWS)
    refs = ReferenceData.from_excel(path)
    with pytest.raises(ReferenceError, match="Unknown nitrate route"):
        refs.get_nitrate_association_dGrxn("some_made_up_route")