"""Reference: every option of nitrate_analysis/calcs.py, with a runnable example.

    python calcs_options.py          (run from the repo root; prints each example)

Copy a block into your own analysis script. Rules behind the options are in
CALCS_RULES.md. Paths come from local_settings.py (DATA_ROOT, LABELS_PATH).
"""

import pandas as pd

from nitrate_analysis.calcs import (
    USABLE_STATES, attach_labels, best_calc, load_calcs, load_labels,
    quality_gate, species_match,
)
from nitrate_analysis.local_settings import DATA_ROOT, LABELS_PATH

pd.set_option("display.width", 200)
pd.set_option("display.max_columns", 30)


def show(title: str, *lines) -> None:
    print(f"\n=== {title}")
    for line in lines:
        print(line)


# ---------------------------------------------------------------------------
# 1. BUILD THE TABLE
# ---------------------------------------------------------------------------
# load_calcs(root, biases=None): one row per calc root (dir containing ion_opt/).
#   biases: tuple of bias dirs to keep, e.g. ("0V", "n1p0V"); None = all.
#   Columns: kind (ads|surf), surf, metal, adsorbate, bias, U_SHE, config,
#            finished, G_eV, F_eV, mu_eV, N, N0, energy_source, out_size, rel
#   finished = ANY file with "finished" in its name in the calc root or in
#              ion_opt/ (finished.txt, sophie_calc_finished.txt, ...). It does
#              NOT check convergence inside `out`.
#   energy_source: cache | out | not parsed (unfinished) | unparseable
calcs = load_calcs(DATA_ROOT, ("0V", "n1p0V"))
show("1a. load_calcs: one row per calculation", f"{len(calcs)} rows", calcs.head(3))

# load_labels(path): {calc rel path: label dict} from ht_labels.json.
labels = load_labels(LABELS_PATH)

# attach_labels(df, labels, require_human_approved=True): adds
#   state       intact | physisorbed | desorbed | broken | None
#   geometry    free text from the label (only some calcs)
#   file_state  e.g. "needs resubmission"
#   label_valid False = no label, stale label (out size changed), or (by default)
#               not human approved. Such rows count as unlabelled.
#   counts_as   ONE species per row (see rules below)
# It never removes calcs. A calc with two pseudonyms becomes two rows (same rel).
#   counts_as rules:
#     intact/physisorbed/desorbed, no pseudonym -> the submitted adsorbate
#     any of those, with pseudonyms             -> the pseudonyms (replace the name)
#     broken, with pseudonyms                   -> the pseudonyms
#     broken without pseudonym, not judged, unlabelled, stale -> None
df = attach_labels(calcs, labels)
# df = attach_labels(calcs, labels, require_human_approved=False)  # use unapproved labels too
show("1b. attach_labels: states and counts_as",
     df["state"].value_counts(dropna=False).to_dict(),
     f"{len(df)} rows from {df['rel'].nunique()} calcs (two-pseudonym calcs are 2 rows)")

# ---------------------------------------------------------------------------
# 2. SELECT: a quality gate AND a species match, combined by you with &
# ---------------------------------------------------------------------------
# quality_gate(df, states=USABLE_STATES, broken_with_pseudonym=True) -> mask
#   finished + has G_eV; adsorbates also need: a counts_as, no file_state flag,
#   and state in `states` (broken passes only if broken_with_pseudonym).
#   Surfaces need no label.
gate = quality_gate(df)
gate_intact_only = quality_gate(df, states=("intact",), broken_with_pseudonym=False)
show("2a. quality_gate", f"default: {gate.sum()} rows pass",
     f"intact only, no broken: {gate_intact_only.sum()} rows pass",
     f"USABLE_STATES = {USABLE_STATES}")

# species_match(df, species, by="counts_as", exclusive=False) -> mask
#   by="counts_as"   strict: what the structure turned out to be
#   by="adsorbate"   loose: what was submitted (dir name), whatever it became
#   exclusive=True   super strict (counts_as only): drop calcs that count as
#                    ANY other species too. Needs the whole table, not a slice.
strict = species_match(df, "NO3")
loose = species_match(df, "NO3", by="adsorbate")
super_strict = species_match(df, "NO3", exclusive=True)
show("2b. species_match for NO3 (before the gate)",
     f"strict {strict.sum()}  loose {loose.sum()}  super strict {super_strict.sum()}")

# Combine them yourself when you want the rows, not just the best one:
usable_no3 = df[gate & strict]
show("2c. gate & species_match", f"{len(usable_no3)} usable NO3 rows",
     usable_no3[["metal", "bias", "config", "state", "G_eV"]].head(3))

# ---------------------------------------------------------------------------
# 3. BEST CALC: lowest G_eV per (metal, bias, species)
# ---------------------------------------------------------------------------
# best_calc(df, metal, bias, species=None, kind="ads", states=USABLE_STATES,
#           by="counts_as", exclusive=False) -> Best(row, status)
#   kind="surf": bare surface (species ignored; includes ASE_pert_from_<bias>).
#   status: "ok" | "no finished calcs" | "no usable calc (broken or unlabelled)"
#   Ties: first rel alphabetically.
def describe(b):
    return (b.status if b.row is None
            else f"{b.row['config']} ({b.row['state']}), G = {b.row['G_eV']:.2f} eV")


show("3. best_calc",
     f"strict       : {describe(best_calc(df, 'Fe', '0V', 'NO3'))}",
     f"loose        : {describe(best_calc(df, 'Fe', '0V', 'NO3', by='adsorbate'))}",
     f"super strict : {describe(best_calc(df, 'Fe', '0V', 'NO3', exclusive=True))}",
     f"intact only  : {describe(best_calc(df, 'Fe', '0V', 'NO3', states=('intact',)))}",
     f"bare surface : {describe(best_calc(df, 'Fe', '0V', kind='surf'))}",
     f"no such calc : {describe(best_calc(df, 'Fe', '0V', 'NOT_A_SPECIES'))}")

# ---------------------------------------------------------------------------
# 4. LOOKING AT WHAT IS THERE (no helper needed, just pandas)
# ---------------------------------------------------------------------------
# Why was a calc not picked? Print its row; compare finished / state /
# label_valid / counts_as / file_state against the gate rules above.
w = df[(df["metal"] == "W") & (df["bias"] == "n1p0V") & (df["adsorbate"] == "NO3")]
show("4a. inspect one metal/bias",
     w[["config", "finished", "state", "counts_as", "label_valid"]].head(5))

# Count calculations (not rows): two-pseudonym calcs are two rows.
show("4b. count calcs, not rows",
     f"rows {len(df)}  calcs {len(df.drop_duplicates('rel'))}")

# Every name a chosen calc counts as: all rows with the same rel.
best = best_calc(df, "Fe", "0V", "NO3").row
if best is not None:
    show("4c. all names of the winning calc",
         df.loc[df["rel"] == best["rel"], "counts_as"].tolist())

# Filter by how a config was made: just use the path/config string.
show("4d. filter by config prefix",
     f"from-scratch (in_*): {df['config'].fillna('').str.startswith('in_').sum()} rows",
     f"from_*: {df['config'].fillna('').str.startswith('from_').sum()} rows")
