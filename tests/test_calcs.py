import json

import pandas as pd
import pytest

from nitrate_analysis.calcs import (
    attach_labels,
    bias_to_volts,
    load_calcs,
    best_calc,
    best_energies,
    candidates,
    gate_reasons,
    load_labels,
    quality_gate,
    species_match,
)

OUT = "ion_opt/jdftx_run/out"


def make_calc(root, rel, *, G=None, out_size=10, finished=True):
    d = root / rel
    run = d / "ion_opt" / "jdftx_run"
    run.mkdir(parents=True)
    (run / "out").write_bytes(b"x" * 10)
    if G is not None:
        (run / ".outfile_cache").mkdir()
        (run / ".outfile_cache" / "data.json").write_text(json.dumps({
            "meta": {"out_size": out_size},
            "data": {"energy": G, "N": 1.0, "N0": 1, "mu": -4.66,
                     "ecomponents": {"F": G - 1}},
        }))
    if finished:
        (d / "ion_opt" / "finished.txt").write_text("Done")
    return d


def label(state, size=10, **extra):
    return {"_out_size": size, "state": state, "human_approved": True, **extra}


@pytest.mark.parametrize("bias, volts", [
    ("0V", 0.0), ("n1p0V", -1.0), ("n0p5V", -0.5), ("1p5V", 1.5), ("No_bias", None),
])
def test_bias_to_volts(bias, volts):
    assert bias_to_volts(bias) == volts


def test_load_calcs_parses_layout_and_cache(tmp_path):
    make_calc(tmp_path, "ads/FeN4-i/NO3/0V/cfg_a", G=-5.0)
    make_calc(tmp_path, "surfs/FeN4-i/0V", G=-3.0)
    make_calc(tmp_path, "surfs/FeN4-i/ASE_pert_from_0V", G=-3.5)
    df = load_calcs(tmp_path)

    ads = df[df["kind"] == "ads"].iloc[0]
    assert (ads["metal"], ads["adsorbate"], ads["bias"], ads["config"]) == ("Fe", "NO3", "0V", "cfg_a")
    assert ads["G_eV"] == -5.0 and ads["F_eV"] == -6.0 and ads["energy_source"] == "cache"
    assert ads["rel"] == "ads/FeN4-i/NO3/0V/cfg_a" and ads["out_size"] == 10

    surfs = df[df["kind"] == "surf"]
    assert set(surfs["bias"]) == {"0V"}  # perturbed re-run counts as 0V
    assert set(surfs["config"].dropna()) == {"ASE_pert_from_0V"}


def test_stale_cache_is_not_trusted(tmp_path):
    make_calc(tmp_path, "surfs/FeN4-i/0V", G=-3.0, out_size=999, finished=False)
    row = load_calcs(tmp_path).iloc[0]
    assert row["energy_source"] == "not parsed"  # unfinished: `out` not parsed
    assert pd.isna(row["G_eV"])


def test_bias_filter(tmp_path):
    make_calc(tmp_path, "surfs/FeN4-i/0V", G=-3.0)
    make_calc(tmp_path, "surfs/FeN4-i/No_bias", G=-3.0)
    assert list(load_calcs(tmp_path, ("0V",))["bias"]) == ["0V"]


def test_load_labels_keeps_ion_opt_and_strips_suffix(tmp_path):
    p = tmp_path / "labels.json"
    p.write_text(json.dumps({"meta": {}, "labels": {
        f"ads/FeN4-i/NO3/0V/a/{OUT}": label("intact"),
        "ads/FeN4-i/NO3/0V/a/sp/jdftx_run/out": label("broken"),
    }}))
    assert load_labels(p) == {"ads/FeN4-i/NO3/0V/a": label("intact")}


def test_attach_labels_ignores_labels_for_a_changed_out(tmp_path):
    make_calc(tmp_path, "ads/FeN4-i/NO3/0V/same", G=-1.0)
    make_calc(tmp_path, "ads/FeN4-i/NO3/0V/rerun", G=-1.0)
    labels = {
        "ads/FeN4-i/NO3/0V/same": label("intact", geometry="bidentate"),
        "ads/FeN4-i/NO3/0V/rerun": label("intact", size=123),  # out changed since
    }
    df = attach_labels(load_calcs(tmp_path), labels).set_index("config")
    assert df.loc["same", "state"] == "intact" and df.loc["same", "geometry"] == "bidentate"
    assert df.loc["same", "counts_as"] == "NO3"
    assert df.loc["rerun", "state"] is None and not df.loc["rerun", "label_valid"]
    assert df.loc["rerun", "counts_as"] is None


def test_unapproved_labels_count_as_unlabelled(tmp_path):
    make_calc(tmp_path, "ads/FeN4-i/NO3/0V/c", G=-1.0)
    labels = {"ads/FeN4-i/NO3/0V/c": label("intact", human_approved=False)}
    df = attach_labels(load_calcs(tmp_path), labels)
    assert df.loc[0, "state"] is None and df.loc[0, "counts_as"] is None
    df = attach_labels(load_calcs(tmp_path), labels, require_human_approved=False)
    assert df.loc[0, "state"] == "intact" and df.loc[0, "counts_as"] == "NO3"


@pytest.mark.parametrize("lab, expected", [
    (label("intact"), ["NO3"]),
    (label("physisorbed"), ["NO3"]),
    (label("desorbed"), ["NO3"]),
    (label("intact", pseudonyms=["O_NO2"]), ["O_NO2"]),  # replaces the name
    (label("broken", pseudonyms=["O_NO2"]), ["O_NO2"]),
    (label("broken"), [None]),
    (label(None), [None]),  # not judged
    (label("intact", size=123), [None]),  # stale
    (None, [None]),  # unlabelled
    (label("intact", pseudonyms=["O_OH_NO2", "O_O_NO2H"]), ["O_OH_NO2", "O_O_NO2H"]),
])
def test_counts_as_rules(tmp_path, lab, expected):
    make_calc(tmp_path, "ads/FeN4-i/NO3/0V/c", G=-1.0)
    labels = {"ads/FeN4-i/NO3/0V/c": lab} if lab else {}
    df = attach_labels(load_calcs(tmp_path), labels)
    assert list(df["counts_as"]) == expected
    assert df["rel"].nunique() == 1  # explode keeps the physical calc's rel


def _gate_df(tmp_path):
    configs = {
        "intact": label("intact"),
        "physisorbed": label("physisorbed"),
        "desorbed": label("desorbed"),
        "broken": label("broken"),
        "broken_pseudo": label("broken", pseudonyms=["O_NO2"]),
        "split": label("intact", pseudonyms=["O_NO2"]),
        "resubmit": label("intact", file_state="needs resubmission"),
        "unlabelled": None,
    }
    labels = {}
    for cfg, lab in configs.items():
        make_calc(tmp_path, f"ads/FeN4-i/NO3/0V/{cfg}", G=-1.0)
        if lab:
            labels[f"ads/FeN4-i/NO3/0V/{cfg}"] = lab
    make_calc(tmp_path, "ads/FeN4-i/NO3/0V/unfinished", G=-1.0, finished=False)
    labels["ads/FeN4-i/NO3/0V/unfinished"] = label("intact")
    make_calc(tmp_path, "surfs/FeN4-i/0V", G=-3.0)  # surfaces need no label
    return attach_labels(load_calcs(tmp_path), labels)


def test_quality_gate_and_species_match(tmp_path):
    df = _gate_df(tmp_path)
    cfgs = lambda m: set(df.loc[m, "config"].fillna("surf"))
    gate = quality_gate(df)
    assert cfgs(gate) == {"intact", "physisorbed", "desorbed", "broken_pseudo", "split", "surf"}
    assert cfgs(quality_gate(df, ("intact",), broken_with_pseudonym=False)) == {
        "intact", "split", "surf"}
    assert cfgs(gate & species_match(df, "NO3")) == {"intact", "physisorbed", "desorbed"}
    assert cfgs(gate & species_match(df, "O_NO2")) == {"broken_pseudo", "split"}
    # loose: everything submitted as NO3 that passes the gate
    assert "split" in cfgs(gate & species_match(df, "NO3", by="adsorbate"))


def test_best_calc_picks_lowest_and_explains_failures(tmp_path):
    for cfg, G in [("in_a", -5.0), ("in_b", -7.0), ("rb_from_broken_X_0", -7.0)]:
        make_calc(tmp_path, f"ads/FeN4-i/NO3/0V/{cfg}", G=G)
    make_calc(tmp_path, "surfs/FeN4-i/0V", G=-3.0)
    make_calc(tmp_path, "surfs/FeN4-i/ASE_pert_from_0V", G=-3.5)
    make_calc(tmp_path, "ads/CoN4-i/NO3/0V/in_a", G=-1.0)
    make_calc(tmp_path, "ads/CoN4-i/NO3/0V/in_u", G=-1.0, finished=False)
    labels = {f"ads/FeN4-i/NO3/0V/{c}": label("intact")
              for c in ("in_a", "in_b", "rb_from_broken_X_0")}
    labels["ads/CoN4-i/NO3/0V/in_a"] = label("broken")
    df = attach_labels(load_calcs(tmp_path), labels)

    best = best_calc(df, "Fe", "0V", "NO3")
    assert best.status == "ok" and best.row["config"] == "in_b"  # tie -> first rel
    assert best_calc(df, "Fe", "0V", kind="surf").row["G_eV"] == -3.5
    assert best_calc(df, "Co", "0V", "NO3").status.startswith("no usable")
    assert best_calc(df, "Zr", "0V", "NO3") == (None, "no finished calcs")
    assert best_calc(df, "Co", "0V", kind="surf").row is None


def test_exclusive_drops_calcs_that_count_as_two_species(tmp_path):
    make_calc(tmp_path, "ads/FeN4-i/NO3/0V/two", G=-9.0)
    make_calc(tmp_path, "ads/FeN4-i/NO3/0V/one", G=-5.0)
    make_calc(tmp_path, "surfs/FeN4-i/0V", G=-3.0)
    labels = {
        "ads/FeN4-i/NO3/0V/two": label("intact", pseudonyms=["NO3", "O_O_NO"]),
        "ads/FeN4-i/NO3/0V/one": label("intact"),
    }
    df = attach_labels(load_calcs(tmp_path), labels)
    assert best_calc(df, "Fe", "0V", "NO3").row["config"] == "two"  # counted in both
    assert best_calc(df, "Fe", "0V", "O_O_NO").row["config"] == "two"
    assert best_calc(df, "Fe", "0V", "NO3", exclusive=True).row["config"] == "one"
    assert best_calc(df, "Fe", "0V", "O_O_NO", exclusive=True).status.startswith("no usable")
    with pytest.raises(ValueError):
        species_match(df, "NO3", by="adsorbate", exclusive=True)


@pytest.mark.parametrize("where, name, expected", [
    ("ion_opt", "finished.txt", True),
    ("", "sophie_calc_finished.txt", True),
    ("", "sophie_finished_sp.txt", True),
    ("ion_opt/jdftx_run", "finished.txt", False),  # subfolders are not looked at
    ("", "adsorbed.txt", False),
])
def test_any_finished_marker_counts(tmp_path, where, name, expected):
    d = make_calc(tmp_path, "surfs/FeN4-i/0V", G=-3.0, finished=False)
    (d / where).mkdir(parents=True, exist_ok=True)
    (d / where / name).write_text("x")
    assert bool(load_calcs(tmp_path).loc[0, "finished"]) is expected


def test_gate_reasons_explain_every_rejection(tmp_path):
    df = _gate_df(tmp_path)
    why = dict(zip(df["config"].fillna("surf") + "/" + df["counts_as"].fillna("-"),
                   gate_reasons(df)))
    assert why["intact/NO3"] is None and why["surf/-"] is None
    assert why["unfinished/NO3"] == "not finished"
    assert why["unlabelled/-"] == "no valid label"
    assert why["broken/-"] == "no counts_as"
    assert why["resubmit/NO3"] == "needs resubmission"
    assert gate_reasons(df, ("intact",), broken_with_pseudonym=False)[
        (df["config"] == "desorbed")].iloc[0] == "state 'desorbed'"
    # the mask and the reasons are the same rule
    assert (quality_gate(df) == gate_reasons(df).isna()).all()


def test_candidates_rank_1_is_best_calc(tmp_path):
    for cfg, G in [("in_a", -5.0), ("in_b", -7.0), ("in_c", -7.0), ("in_d", -9.0)]:
        make_calc(tmp_path, f"ads/FeN4-i/NO3/0V/{cfg}", G=G)
    labels = {f"ads/FeN4-i/NO3/0V/{c}": label("intact") for c in ("in_a", "in_b", "in_c")}
    labels["ads/FeN4-i/NO3/0V/in_d"] = label("broken")  # lowest energy but cannot win
    df = attach_labels(load_calcs(tmp_path), labels)

    cand = candidates(df, "Fe", "0V", "NO3")
    assert cand.loc[cand["rank"] == 1, "rel"].item() == best_calc(df, "Fe", "0V", "NO3").row["rel"]
    assert cand["config"].tolist() == ["in_b", "in_c", "in_a", "in_d"]  # passers first
    assert cand["dG_vs_best"].round(6).tolist() == [0.0, 0.0, 2.0, -2.0]
    assert cand["why_not"].tolist() == [None, None, None, "no counts_as"]
    assert cand["rank"].isna().tolist() == [False, False, False, True]


def test_candidates_surface_and_no_match(tmp_path):
    make_calc(tmp_path, "surfs/FeN4-i/0V", G=-3.0)
    make_calc(tmp_path, "surfs/FeN4-i/ASE_pert_from_0V", G=-3.5)
    df = attach_labels(load_calcs(tmp_path), {})
    surf = candidates(df, "Fe", "0V", kind="surf")
    assert surf["G_eV"].tolist() == [-3.5, -3.0] and surf["rank"].tolist() == [1, 2]
    assert candidates(df, "Fe", "0V", "NO3").empty
    with pytest.raises(ValueError):
        candidates(df, "Fe", "0V")


def test_best_energies_returns_energies_and_inspectable_picks(tmp_path):
    make_calc(tmp_path, "surfs/FeN4-i/0V", G=-3.0)
    make_calc(tmp_path, "ads/FeN4-i/NO3/0V/in_a", G=-5.0)
    make_calc(tmp_path, "ads/FeN4-i/NO2/0V/in_a", G=-4.0, finished=False)
    df = attach_labels(load_calcs(tmp_path), {"ads/FeN4-i/NO3/0V/in_a": label("intact")})
    energies, picks = best_energies(df, "Fe", "0V", {"*": None, "*NO3": "NO3", "*NO2": "NO2", "*X": "X"})
    assert energies == {"*": -3.0, "*NO3": -5.0}  # keys without a usable calc are left out
    assert picks["*NO3"].row["config"] == "in_a"
    assert picks["*NO2"].status == "no finished calcs" and picks["*X"].status == "no finished calcs"
    # kwargs reach best_calc
    assert "*NO3" not in best_energies(df, "Fe", "0V", {"*NO3": "NO3"}, states=("desorbed",))[0]
