# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Purpose

Post-processing and thermodynamic analysis of JDFTx DFT calculations for the
nitrate reduction reaction (NO3RR) on M-N-C (metal–nitrogen–carbon) catalysts.
This is research code, not a distributed package.

## Commands

```bash
# Run the full test suite (run from the repo root so `nitrate_analysis` imports)
python -m pytest -q

# Single file / test / pattern
python -m pytest tests/test_fed.py -q
python -m pytest tests/test_valence.py::test_known_formulas
python -m pytest -k "pcet"
```

There is no packaging (`pyproject.toml`/`setup.py`), lint, or formatter config.
Dependencies are implicit: `pandas`, `openpyxl` (Excel engine), `matplotlib`, `pytest`, and
`pymatgen` (only to parse a JDFTx `out` file when the cache is missing/stale).
Imports resolve only because pytest/scripts are run from the repo root.

## Working style

* **Core code lives in `nitrate_analysis/`; scripts at the repo root are
  standalone and must not import each other.** A script (`get_Fe_FED.py`,
  `compare_no3_adsorption.py`, `calcs_options.py`) imports only from
  `nitrate_analysis/`. Anything two scripts both need moves into a core module.
* Priority is the core modules. The root scripts are throwaway/analysis
  drivers and are not kept in sync with API changes unless asked.

## Architecture

Six cooperating modules under `nitrate_analysis/`, each with one job. Data
flows `calcs` (energies) + `references` (constants) + `valence` + `network` (graph) → `fed`. `fed`
knows nothing about files or the workbook; the calling script wires them.

| module | responsibility |
|---|---|
| `calcs.py` | Read the local copy of the calculation tree into one DataFrame (`load_calcs`), attach the collaborator's hand-checked labels (`load_labels`, `attach_labels`), and pick calcs: `quality_gate`, `species_match`, `best_calc(df, metal, bias, species, kind=...) → Best(row, status)` = lowest `G_eV` among gated rows; `best_energies(df, metal, bias, {key: species or None})` → `(energies, picks)` runs it for many species at once. `gate_reasons` is the one definition of the gate (why a row fails; `quality_gate` = no reason), and `candidates(...)` lists every calc `best_calc` considered with its rank and `why_not`. Knows nothing about reactions or reference energies. Rules and open questions: `CALCS_RULES.md`. Usage examples for every option: `calcs_options.py`. |
| `references.py` | `ReferenceData.from_excel(path)` loads the cached `mu_*` / `dG_rxn_*` values from `References.xlsx` (sheet `ref_summary`, columns `name`, `value_eV`). `get(name)` raises `ReferenceError` (lists available names) for an unknown name. Also exposes two *derived* names, `mu_e0` (= `MU_E0_SHE`, −4.66 eV) and `mu_H+` (= ½·mu_H2 − mu_e0, pH 0). No DFT here. |
| `valence.py` | `valence_electrons("NO3-") → 24`: valence electron count from a formula (optional trailing `+`/`-`), using GBRV counts in `ZVAL` (H, N, O so far). Raises `ValueError` for unknown elements/unparseable formulas. Add metals to `ZVAL` from the pseudopotential files when needed (GBRV Fe is 16, not 8). |
| `network.py` | The reaction network as data: `NO3RR_GRAPH` (edge → per-step `n_h2o`/`n_h`, optional `"ref"`). Currently two routes NO3- → \*NO2 (then \*NO2 → NO2- desorption): intact (\*NO3 → \*NO3H → \*NO2) and split (\*O_NO2 → \*NO2). Adsorbed state names equal the calc species names after dropping `*`. |
| `fed.py` | Free-energy diagrams, no data of its own. `find_paths`/`find_path` walk a graph; `graph_states(graph)` lists the states whose energy must come from calcs; `get_FEDs(metal, start, end, bias, *, graph, energies, constants)` prices every route (`get_FED` = first route only) with `ion_adsorption_dG` (aqueous ↔ `*X`) or `pcet_dG` (`*X → *Y`, optional `"ref"` state the dG is measured against). |
| `plotting.py` | `plot_FED(routes, *, zero="initial"\|"final", route_labels, ax, title) → Axes` draws `get_FEDs` output (bars + dashed connectors, `n/a` gap for a `None` energy). `zero` re-zeros each route at its start or end state (`shift_route`). Knows nothing about files/calcs; the caller saves or shows. |

`script.py` is an unused stub (imports modules that no longer exist).

### fed.py conventions

* `graph`, `energies` (state → G in eV, `"*"` = bare slab) and `constants` are
  **required**, no defaults, so fake numbers can't leak into a real result.
  `constants` is a `ReferenceData` *or* a plain dict, using the workbook's
  names: `mu_e0`, `mu_H+`, `mu_H2O`, `mu_HNO3`, `mu_NO3_H`, `mu_NO2_H` (chemical potential of the ion + H+ pair), `dG_rxn_redox`.
* **Two kinds of "missing":** a missing *constant* is a setup bug and raises
  (`ReferenceError` from `ReferenceData`, `KeyError` from a dict). A missing
  *energy* gives `None` for that state (printed `--` in a diagram); it only
  blocks states whose edge actually references it.
* Bias strings look like `"0V"`, `"-1.0V"`, `f"{U_SHE}V"`. The calc-tree
  directory name `n1p0V` is *not* parseable by `fed`; use `calcs`' `U_SHE`.
* `mu_e = mu_e0 − U`. `mu_e0` is hardcoded: the JDFTx fixed-potential runs
  report −4.660 eV at 0V (spread ~1e-5) and exactly 1 eV higher at −1V.
* Graph states are case-insensitive but canonicalised to uppercase, so
  `energies` keys must be uppercase (species names in the calc tree are).
* Routes that end in the same state must agree on its energy (intermediates
  cancel); a mismatch means a step's stoichiometry in the graph is wrong.
* The graph ends at `NO2-`: both routes reach `*NO2`, then desorb via `*NO2 → NO2-` (`desorption=True`).
* `ion_adsorption_dG(ion, surface_state, bias, ..., desorption=False)` is
  `G_surf − G_* − mu_e·(−N_ion) − mu_ion_H+ + mu_H+` for `ion(aq) + * → *X + e-`,
  with `N_ion` from `valence_electrons(ion)` and `mu_<ion>_H` the workbook name
  `mu_<formula>_H`, charge dropped (`mu_NO3_H`, `mu_NO2_H`). `desorption=True` returns the
  negative (`*X + e- → ion(aq) + *`). `fed` picks adsorption vs desorption from
  which end of an edge is aqueous. `nitrate_adsorption_dG(to_state, ...)` is a
  thin wrapper for `ion="NO3-"` (no `dG_rxn_redox`). The old formula, which adds `dG_rxn_redox`, is kept as
  `nitrate_adsorption_dG_cycle` for reference only; nothing uses it.

### Tests

`tests/` mirrors the modules. `tests/example_data.py` holds the **fake**
energies/constants for `fed` tests (`test_fed.py` binds them with `partial`).
Fake data belongs only under `tests/`.

### Known issue: workbook `value_eV` may be Hartree

Historically the `ref_summary` `mu_*` magnitudes looked like Hartree (e.g.
`mu_H2 ≈ −31.72`) while the calc energies are eV. Not re-verified after the
workbook changes. Any diagram mixing the two scales is wrong; check before
trusting absolute energies, and fix in one place (`references.py` load).

## Config / secrets

`nitrate_analysis/local_settings.py` (git-ignored) holds machine-local absolute
paths: `WORKBOOK_PATH` (`References.xlsx`), `DATA_ROOT` (calc tree copy),
`LABELS_PATH` (`ht_labels.json`). Never put absolute paths in tracked modules.
