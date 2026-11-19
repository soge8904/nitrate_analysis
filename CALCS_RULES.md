# calcs.py: rules and conventions

Working spec for `nitrate_analysis/calcs.py` (the data layer). Nothing here is
implemented yet unless it says "current". Items marked **⚠ CONFIRM** are my
reading or a guess and need your check. Items marked **TODO** are decisions
still open. Edit freely.

Data facts below come from `ht_labels.json` (exported 2026-11-06, 6754 labels).

---

## 1. Role of calcs.py

- The data layer: filesystem + labels in, a table out. It answers "what
  calculations exist, and which can I trust?"
- It never prices anything. No reactions, no reference energies, no chemistry.
  Anything physical is passed in as an argument.
- Layers: **load** (calcs.py) -> **select** (calcs.py) -> **price**
  (fed.py, references.py) -> **analyse / plot** (standalone scripts).
- Analysis scripts should contain no filtering code of their own, only
  choices (which metals, biases, species, states).

## 2. Directory layout and row identity (current)

One row per calculation root, i.e. a directory containing `ion_opt/`.

| kind | path | notes |
|---|---|---|
| `ads` | `ads/<surf>/<adsorbate>/<bias>/<config>/` | 5 path parts |
| `surf` | `surfs/<surf>/<bias_dir>/` | 3 path parts |

- `metal` = `surf.split("N4")[0]` (e.g. `CoN4-i` -> `Co`).
- `bias` is `0V`, `n1p0V` (= -1.0 V vs SHE), `n0p5V`, or `No_bias`.
  `No_bias` is a fixed-charge run: `U_SHE` is `None`. Only `0V` and `n1p0V`
  are fixed-potential runs used for comparison.
- `rel` = the path relative to the data root. It is the key into the labels
  file.
- `DATA_ROOT` currently points at the NO3 subset
  (`data/no3_subset`), not the full tree. Checks that need other adsorbates
  must use the labels file (which lists everything) or the full tree.

## 3. Energies (current)

- `finished` = `ion_opt/finished.txt` exists.
- `G_eV` is JDFTx's grand free energy `F - mu*N` for fixed-potential runs.
  `mu_eV`, `N` are the electron chemical potential and electron count.
- Source: `ion_opt/jdftx_run/.outfile_cache/data.json` if
  `meta.out_size` == current size of `out`; otherwise parse `out` with
  pymatgen's `JDFTXOutfile`. Unfinished calcs are not parsed (slow, ~0.2 s).
- `energy_source` records which path was used: `cache`, `out`,
  `not parsed`, `unparseable`.

## 4. Config-name prefixes

The directory name under `<bias>/` says how the starting geometry was made.

| prefix | meaning | approx. count (ads) |
|---|---|---|
| `in_` | built from scratch | 4739 |
| `from_<...>` | initial geometry = the **final** geometry of another calc, with a modification. Example: `O_NH/0V/from_o_n_gau_gau_anti` = final geometry of `O_N/0V/in_gau_gau_anti` with an H added to N | 1802 |
| `ASE_pert...` | imaginary modes of the original calc were perturbed and the structure re-optimised. For bare surfaces the directory is `ASE_pert_from_<bias>` | 51 |
| `ASE_from...` | a **normal** calc that used the ASE optimiser. Not a perturbation | ~54 |
| `rb_from_broken_<adsorbate>_<n>` | the broken/changed-species case, see section 6. **Not** the same as `from_` | 93 |
| `man1` | ignore | 11 |

Other things seen in names:

- Suffix `_resub` (376 names). **⚠ CONFIRM**: I read it as "resubmitted run".
  Does it affect how the calc should be treated?
- `from_` names do not encode the source consistently
  (`from_o_n_gau_gau_anti`, `from_OH_OH_NO2H_in_bi_gau_sym_n_coord_a`,
  `from_n1p0V_bi_gau_sym_n_coord`). Do not parse the source out of the name.
- `from_` chains exist (a `from_` calc can itself be the source of another).
- **⚠ CONFIRM**: `rb` stands for what? (I do not know; do not guess in code.)

Surfaces (current): a directory named `ASE_pert_from_<bias>` gets
`bias=<bias>` and `config=<dirname>`, so it competes with the plain surface
for "lowest energy at that bias". Adsorbate configs are stored as the plain
`config` string with no special handling.

**TODO**: add an `origin` column derived from `config`, with values
`scratch | from | pert | ase | rb | other`, so scripts can say "only
from-scratch calcs" without string matching. Needs your approval of the table
above first. `man1` handling: drop the rows, or `other`?

## 5. Labels (ht_labels.json)

- Labels are the collaborator's hand-checked judgments and **supersede** the
  `intact.txt` / `adsorbed.txt` / `broken.txt` marker files in the calc roots.
  Do not fall back on the marker files.
- A label is valid only if its `_out_size` equals the current `out` size. A
  rerun changes the size; a stale label is treated as **unlabelled**.
- Calcs with no label are left out of the collaborator's analysis. Never
  usable (for adsorbates).
- Fields: `state`, `geometry` (only some), `pseudonyms` (only some),
  `file_state`, `human_approved`.
- `state` counts: intact 4902, physisorbed 716, desorbed 609, broken 500,
  not judged 27.
- `file_state = "needs resubmission"` -> not usable.
- `human_approved` is true for 6733 of 6754. **⚠ CONFIRM**: ignore it, or
  treat the 21 unapproved labels specially? (Current code ignores it.)
- The collaborator's own usable states are intact / physisorbed / desorbed
  (the "usable states" in the file's meta). A stricter script of theirs uses
  only intact / physisorbed.
- `desorbed` is usable: a desorbed config is still a valid energy for the
  slab + species composition.

## 6. Pseudonyms, `broken`, and `counts_as`

`adsorbate` = the directory name = **what was submitted**.
`counts_as` = **what the structure turned out to be** (from the label).

### Rules (decided)

| label | `counts_as` | usable? |
|---|---|---|
| valid, intact / physisorbed / desorbed, no pseudonym | the adsorbate | yes |
| valid, intact / physisorbed / desorbed, with pseudonyms | the pseudonyms | yes, as the pseudonyms |
| valid, **broken**, with pseudonyms | the pseudonyms | **yes** (bonds broke, fragments are still adsorbed) |
| valid, **broken**, no pseudonym | `None` | **no** |
| valid label, no state ("not judged") | `None` | no |
| missing or stale label | `None` | no |
| `file_state` = needs resubmission | (unchanged) | no |

- `broken` does **not** mean unusable. It means a bond in the calc broke.
  If every fragment is still adsorbed it is usable and carries a pseudonym. No
  pseudonym means not usable.
- Pseudonyms **replace** the default name. An NO3 calc with pseudonym
  `['O_NO2']` counts as `O_NO2` and not as `NO3`.
- Current code gets this wrong in one place: for any valid label without
  pseudonyms, `counts_as` defaults to `[adsorbate]`, so a broken NO3 with no
  pseudonym claims to count as NO3 (432 broken + 27 not-judged rows are
  affected). Must be fixed.
- Pseudonym lists have 1 or 2 names (never more). About 95 labels have two.
  Two names = the structure is ambiguous between two intermediates, e.g.
  `('O_OH_NO2', 'O_O_NO2H')` (H on a different O). One label is
  `('NO3', 'O_O_NO')`.
- 96 pseudonym entries name the calc's own adsorbate (no species change).
  Those are not "turned into something else".
- Broken-with-pseudonym: 68 labels. Examples: `NO2` calcs broken with
  pseudonym `O_NO` (TaN4); `OH_OH_NO2H` broken with `OH_H2O_NO2` (ReN4).

### `rb_from_broken` copies

When a calc turns into another species, its directory is copied under the new
species' directory with the name `rb_from_broken_<source adsorbate>_<n>`.
Example: `ads/ZrN4-i/OH_H2O_NO2/0V/rb_from_broken_OH_OH_NO2H_0`.

- The copy carries the pseudonym label too, so it appears under the species
  it turned into.
- **TODO (check requested)**: for every calc with a pseudonym naming a
  different species, does a matching `rb_from_broken_...` copy exist under the
  new species? Preliminary, using only identical `out` file size under the
  target species (any name, not only `rb_from_broken`): 87 of 168
  (calc, new-species) pairs have a copy; 81 do not (54 broken, 27 intact).
  A re-run copy would have a different size, so 81 is an upper bound on
  "missing". Needs re-running restricted to `rb_from_broken` names, ideally on
  the full tree.
- **⚠ CONFIRM**: the original and its `rb_` copy are the same physical
  calculation. When counting calculations, de-duplicate them. For "lowest
  energy" it should not matter (same energy), but the original has the
  `adsorbate` of the old species. Is that right?

## 7. Table shape: one row per (calculation, species)

- `load_calcs` stays one row per calculation (it only reads disk).
- `attach_labels` makes `counts_as` a single string per row and
  **explodes** calcs with two pseudonyms into two rows (same `rel`, same
  `geometry`, same energy; different `counts_as`).
  - Rows with `counts_as = None` stay as a single row.
- Consequences:
  - `rel` is no longer unique. A row's identity is `(rel, counts_as)`;
    "the physical calc" is `rel`.
  - Anything that **counts calculations** must `drop_duplicates("rel")`
    first. E.g. `compare()`'s `n_finished` would otherwise double-count the
    `('NO3', 'O_O_NO')` calc.
  - Anything that takes a lowest energy **per species** is safe (the
    duplicates sit in different species groups).
  - To see every name a winner counts as: look up all rows with the same
    `rel`.
  - Keep `adsorbate` and `counts_as` as separate columns. Never merge them.

## 8. Selecting calcs: two independent filters

Split `usable()` into two pieces that the caller composes:

1. **Quality gate**: finished, has an energy, valid label, not flagged for
   resubmission, and the section-6 rule (state is intact / physisorbed /
   desorbed, or broken **with** a pseudonym). Allowed states are an argument,
   e.g. `("intact", "physisorbed")` to drop desorbed.
2. **Species match**, chosen by the caller:
   - strict: `counts_as == X` ("NO3 stays NO3")
   - loose: `adsorbate == X` ("everything run as NO3, whatever it became")
   - all intermediates: each `counts_as == X` in turn

| analysis wants | species filter | notes |
|---|---|---|
| NO3 that stayed NO3 | `counts_as == "NO3"` | excludes calcs that became O_NO2 etc. |
| everything run as NO3, even if it changed | `adsorbate == "NO3"` | **TODO**: do unlabelled and broken-no-pseudonym calcs belong here? The current gate requires a valid label, so they would be excluded. |
| every state of intermediate X | `counts_as == X` | a calc can appear under two X's (section 7) |

Surfaces need only finished + energy (no labels).

## 9. Selection helper (new)

`best_calc(df, metal, bias, species, kind, states, by="counts_as")`:

- lowest `G_eV` among rows passing the gate and the species filter;
- returns the row, or `None`;
- should say **why** nothing was found, using the reasons `compare()` already
  reports: no bare surface / no finished calcs / no usable calc (broken or
  unlabelled);
- the bare surface is the lowest `G_eV` among usable `surf` rows at that
  metal and bias, including `ASE_pert_from_<bias>` re-optimisations;
- **TODO**: return type (row, energy, or a small result object carrying the
  reason).
- **TODO**: tie-breaking when two rows have the same energy (e.g. a calc and
  its `rb_` copy).

## 10. Changes this implies (not yet made)

- `attach_labels`: new `counts_as` rule (section 6), explode (section 7).
- `usable`: split into gate + species match (section 8); remove the
  `.map(lambda ...)` list check.
- Add an `origin` column (section 4).
- New `best_calc` (section 9); `compare_no3_adsorption.compare()` calls it
  instead of its own masks, and counts with `drop_duplicates("rel")`.
- `tests/test_calcs.py`: assumes list-valued `counts_as` and one row per calc;
  update it and add tests for each row of the section-6 table.

## 11. Open questions (collected)

1. `rb` meaning and whether `rb_` copies should be de-duplicated (6).
2. Does `_resub` change how a calc is treated? (4)
3. Ignore `human_approved`, or treat the 21 unapproved labels specially? (5)
4. Loose mode: include unlabelled / broken-no-pseudonym calcs? (8)
5. Is the full calc tree available locally for the `rb_from_broken` check? (6)
6. `man1` rows: drop or `other`? (4)
7. Should `broken`-with-pseudonym be selectable by default, or only when
   asked for? (My earlier suggestion was opt-in; your rule now says usable.)
