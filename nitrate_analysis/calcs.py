"""Read a local copy of the NO3_ads_gen calculation tree into one DataFrame.

One row per calculation root (a directory containing ``ion_opt/``):

* ``ads/<surf>/<adsorbate>/<bias>/<config>/``  -> kind ``"ads"``
* ``surfs/<surf>/<bias_dir>/``                 -> kind ``"surf"``

``<bias_dir>`` is normally the bias (``0V``, ``n1p0V``, ``No_bias``), but a few
surfaces also have perturbed re-optimisations named ``ASE_pert_from_<bias>``;
those get ``bias=<bias>`` and ``config=<bias_dir>`` so they compete with the
plain surface for "lowest energy" at that bias.

Energies come from ``ion_opt/jdftx_run/.outfile_cache/data.json`` when the
cache is valid (``meta.out_size`` == size of ``out``); otherwise ``out`` is
parsed with pymatgen's ``JDFTXOutfile`` (the same parser that wrote the cache).
``G_eV`` is JDFTx's grand free energy ``F - mu*N`` for fixed-potential runs;
``mu_eV`` and ``N`` are the electron chemical potential and electron count.

Adsorbate calcs are judged by the collaborator's hand-checked labels
(``ht_labels.json``, see :func:`load_labels`), not by the ``intact.txt`` /
``adsorbed.txt`` marker files in the calc roots: the labels supersede them.

Table shape: ``load_calcs`` gives one row per calculation; ``attach_labels``
gives one row per (calculation, species) -- a calc with two pseudonyms becomes
two rows sharing ``rel``. Count calculations with ``drop_duplicates("rel")``.
Selection = :func:`quality_gate` + a species match, composed by the caller (see
:func:`best_calc`). Rules: CALCS_RULES.md.

Nothing here knows about reactions or reference energies.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import NamedTuple

import pandas as pd

COLUMNS = [
    "kind", "surf", "metal", "adsorbate", "bias", "U_SHE", "config",
    "finished", "G_eV", "F_eV", "mu_eV", "N", "N0", "energy_source",
    "out_size", "rel",
]
LABEL_COLUMNS = ["state", "geometry", "counts_as", "file_state", "label_valid"]

# States the collaborator's pathway analysis treats as usable (ht_labels.json
# meta.usable_states). "desorbed" is in: a desorbed NO3 config is still a valid
# energy for the slab + NO3 composition.
USABLE_STATES = ("intact", "physisorbed", "desorbed")
_PERT_PREFIX = "ASE_pert_from_"
_OUT = "ion_opt/jdftx_run/out"


def bias_to_volts(bias: str) -> float | None:
    """``"0V"`` -> 0.0, ``"n1p0V"`` -> -1.0, ``"n0p5V"`` -> -0.5,
    ``"No_bias"`` -> None (fixed-charge run, no applied potential)."""
    if bias == "No_bias":
        return None
    s = bias.removesuffix("V")
    sign = -1.0 if s.startswith("n") else 1.0
    return sign * float(s.lstrip("n").replace("p", "."))


def _energies(run: Path, parse_out: bool = True) -> dict:
    """G/F/mu/N from the cache if valid, else (if ``parse_out``) by parsing
    ``out``; {} if there is no ``out``."""
    out = run / "out"
    cache = run / ".outfile_cache" / "data.json"
    if not out.exists():
        return {}
    size = out.stat().st_size
    try:
        c = json.loads(cache.read_text())
        if c["meta"]["out_size"] == size:
            d = c["data"]
            return {
                "G_eV": d["energy"], "F_eV": d.get("ecomponents", {}).get("F"),
                "mu_eV": d["mu"], "N": d["N"], "N0": d["N0"],
                "energy_source": "cache", "out_size": size,
            }
    except (FileNotFoundError, KeyError, ValueError):
        pass
    if not parse_out:
        return {"energy_source": "not parsed", "out_size": size}
    try:
        from pymatgen.io.jdftx.outputs import JDFTXOutfile

        o = JDFTXOutfile.from_file(out)
        return {
            "G_eV": o.e, "F_eV": None, "mu_eV": o.mu, "N": o.total_electrons,
            "N0": None, "energy_source": "out", "out_size": size,
        }
    except Exception:  # unparseable/truncated out: leave energies empty
        return {"energy_source": "unparseable", "out_size": size}


def _is_finished(d: Path) -> bool:
    """A calc counts as finished (= converged) if its root or its ``ion_opt/``
    holds ANY file with "finished" in the name: ``finished.txt``,
    ``sophie_calc_finished.txt``, ``sophie_finished_sp.txt``,
    ``finished.txt.pre_unfreeze``. Subfolders (``sp/``, ``vib/``) are not
    looked at. The ``out`` file is not checked for convergence. To be stricter,
    change the name patterns here."""
    return any(any(d.glob("*finished*")) for d in (d, d / "ion_opt"))


def _row(root: Path, d: Path, biases: tuple[str, ...] | None) -> dict | None:
    rel = d.relative_to(root).parts
    if rel[0] == "ads" and len(rel) == 5:
        _, surf, ads, bias, config = rel
        kind = "ads"
    elif rel[0] == "surfs" and len(rel) == 3:
        _, surf, bias_dir = rel
        kind, ads = "surf", None
        if bias_dir.startswith(_PERT_PREFIX):
            bias, config = bias_dir.removeprefix(_PERT_PREFIX), bias_dir
        else:
            bias, config = bias_dir, None
    else:
        return None
    if biases is not None and bias not in biases:
        return None
    row = {
        "kind": kind, "surf": surf, "metal": surf.split("N4")[0],
        "adsorbate": ads, "bias": bias, "U_SHE": bias_to_volts(bias),
        "config": config,
        "finished": _is_finished(d),
        "rel": "/".join(rel),
    }
    # parsing `out` is slow (~0.2 s); only worth it for finished calcs
    row.update(_energies(d / "ion_opt" / "jdftx_run", parse_out=row["finished"]))
    return row


def load_calcs(root: str | Path, biases: tuple[str, ...] | None = None) -> pd.DataFrame:
    """All calc roots under ``root`` (the local PBED30_GBRV_20Ha copy),
    optionally only those at ``biases`` (e.g. ``("0V", "n1p0V")``)."""
    root = Path(root)
    rows = []
    for dirpath, dirnames, _ in os.walk(root):
        if "ion_opt" in dirnames:
            r = _row(root, Path(dirpath), biases)
            if r is not None:
                rows.append(r)
            dirnames.clear()  # a calc root's subdirs are never calc roots
    return pd.DataFrame(rows, columns=COLUMNS)


def load_labels(path: str | Path) -> dict[str, dict]:
    """``{calc root rel path: label}`` from ``ht_labels.json``, ion_opt labels
    only (a few labels sit on sp/ or vib/ out files; the energies here are
    ion_opt's)."""
    labels = json.loads(Path(path).read_text())["labels"]
    suffix = "/" + _OUT
    return {k.removesuffix(suffix): v for k, v in labels.items() if k.endswith(suffix)}


def _counts_as(label: dict, adsorbate: str | None) -> list[str | None]:
    """Species a calc counts as (CALCS_RULES.md section 6)."""
    state = label.get("state")
    pseudo = list(label.get("pseudonyms") or [])
    if state in USABLE_STATES:
        return pseudo or [adsorbate]
    if state == "broken" and pseudo:
        return pseudo  # fragments are still adsorbed
    return [None]  # broken without pseudonym, not judged, unlabelled, stale


def attach_labels(df: pd.DataFrame, labels: dict[str, dict],
                  require_human_approved: bool = True) -> pd.DataFrame:
    """Copy of ``df`` with the label columns, one row per
    (calculation, species):

    * ``state``: intact | physisorbed | desorbed | broken, or None (unlabelled)
    * ``geometry``: e.g. monodentate / bidentate (only some adsorbates)
    * ``counts_as``: the single species this row counts as. Pseudonyms replace
      the submitted ``adsorbate``; broken with pseudonym counts as the
      pseudonym; broken without one, not judged, unlabelled -> None. A calc
      with two pseudonyms is exploded into two rows (same ``rel``).
    * ``file_state``: e.g. "needs resubmission"
    * ``label_valid``: the label is usable: its ``_out_size`` equals the current
      ``out`` size (else the run changed since labelling) AND, when
      ``require_human_approved``, ``human_approved`` is true. An invalid label
      is treated as unlabelled.
    """
    df = df.copy()
    lab = [labels.get(r) for r in df["rel"]]
    valid = [l is not None and l.get("_out_size") == s for l, s in zip(lab, df["out_size"])]
    # HUMAN APPROVAL: labels without human_approved == True count as unlabelled.
    # Pass require_human_approved=False to use them.
    if require_human_approved:
        valid = [v and l.get("human_approved") is True for l, v in zip(lab, valid)]
    use = [l if v else {} for l, v in zip(lab, valid)]
    df["state"] = [l.get("state") for l in use]
    df["geometry"] = [l.get("geometry") for l in use]
    df["counts_as"] = [_counts_as(l, a) for l, a in zip(use, df["adsorbate"])]
    df["file_state"] = [l.get("file_state") for l in use]
    df["label_valid"] = valid
    df = df.explode("counts_as", ignore_index=True)
    df["counts_as"] = df["counts_as"].astype(object).where(df["counts_as"].notna(), None)
    return df


def gate_reasons(df: pd.DataFrame, states: tuple[str, ...] = USABLE_STATES,
                 broken_with_pseudonym: bool = True) -> pd.Series:
    """Why each row fails the quality gate, or None if it passes. The first
    failing rule wins, in this order: ``not finished``, ``no energy``, then
    for adsorbates ``no valid label``, ``state '<state>'`` (not in ``states``),
    ``needs resubmission``, ``no counts_as``. This is the single definition of
    the gate; :func:`quality_gate` is just "no reason"."""
    allowed = set(states) | ({"broken"} if broken_with_pseudonym else set())
    reasons = []
    for r in df.itertuples(index=False):
        if not r.finished:
            reasons.append("not finished")
        elif pd.isna(r.G_eV):
            reasons.append("no energy")
        elif r.kind != "ads":
            reasons.append(None)
        elif r.state is None or pd.isna(r.state):
            reasons.append("no valid label")
        elif r.state not in allowed:
            reasons.append(f"state {r.state!r}")
        elif not pd.isna(r.file_state):
            reasons.append(r.file_state)
        elif r.counts_as is None or pd.isna(r.counts_as):
            reasons.append("no counts_as")
        else:
            reasons.append(None)
    return pd.Series(reasons, index=df.index, dtype=object)


def quality_gate(df: pd.DataFrame, states: tuple[str, ...] = USABLE_STATES,
                 broken_with_pseudonym: bool = True) -> pd.Series:
    """Mask of rows trustworthy enough to price (species-independent).

    Every row: finished and has an energy. Surfaces need nothing else.
    Adsorbate rows (needs :func:`attach_labels`): valid label, not flagged for
    resubmission, ``counts_as`` set, and state in ``states`` -- or ``broken``
    (which only has ``counts_as`` if it carries a pseudonym) when
    ``broken_with_pseudonym``. See :func:`gate_reasons` for the per-row reason.
    """
    return gate_reasons(df, states, broken_with_pseudonym).isna()


def species_match(df: pd.DataFrame, species: str, by: str = "counts_as",
                  exclusive: bool = False) -> pd.Series:
    """Mask for the species filter. ``by="counts_as"``: what the structure
    turned out to be (strict). ``by="adsorbate"``: what was submitted (loose).

    ``exclusive`` (super strict, only with ``by="counts_as"``): also require
    that the calc counts as *only* ``species``. A calc with two pseudonyms
    (two rows sharing ``rel``) is dropped. ``df`` must be the whole table, not
    a pre-filtered slice, or the other row of a two-name calc is not seen.
    """
    if by not in ("counts_as", "adsorbate"):
        raise ValueError(f"by must be 'counts_as' or 'adsorbate', not {by!r}")
    hit = df[by] == species
    if exclusive:
        if by != "counts_as":
            raise ValueError("exclusive only applies to by='counts_as'")
        hit &= df.groupby("rel")["counts_as"].transform("nunique") == 1
    return hit


class Best(NamedTuple):
    """``row`` is the lowest-``G_eV`` row (or None); ``status`` is "ok" or why
    there is none: "no finished calcs" | "no usable calc (broken or
    unlabelled)"."""
    row: pd.Series | None
    status: str


def best_calc(df: pd.DataFrame, metal: str, bias: str, species: str | None = None,
              kind: str = "ads", states: tuple[str, ...] = USABLE_STATES,
              by: str = "counts_as", exclusive: bool = False) -> Best:
    """Lowest ``G_eV`` among gated rows for ``metal`` at ``bias``.

    ``kind="surf"`` is the bare surface (``species`` ignored; includes
    ``ASE_pert_from_<bias>`` re-optimisations). ``kind="ads"`` needs ``species``
    and matches it by ``by`` (and ``exclusive``, see :func:`species_match`). Ties go to the first ``rel`` alphabetically.
    """
    base = (df["kind"] == kind) & (df["metal"] == metal) & (df["bias"] == bias)
    here = base
    if kind == "ads":
        if species is None:
            raise ValueError("species is required for kind='ads'")
        # "no finished calcs" looks at everything run as, or counting as, species
        related = base & ((df["adsorbate"] == species) | species_match(df, species, by))
        if not df.loc[related, "finished"].astype(bool).any():
            return Best(None, "no finished calcs")
        here = base & species_match(df, species, by, exclusive)
    elif not df.loc[base, "finished"].astype(bool).any():
        return Best(None, "no finished calcs")
    cand = df[here & quality_gate(df, states)]
    if cand.empty:
        return Best(None, "no usable calc (broken or unlabelled)")
    return Best(cand.sort_values(["G_eV", "rel"]).iloc[0], "ok")


def best_energies(df: pd.DataFrame, metal: str, bias: str,
                  species: dict[str, str | None], **kwargs
                  ) -> tuple[dict[str, float], dict[str, Best]]:
    """:func:`best_calc` for many species at once. ``species`` maps a key of
    your choosing to a species name (``None`` = the bare surface), e.g.
    ``{"*": None, "*NO3": "NO3"}``; ``kwargs`` (``states``, ``by``,
    ``exclusive``) go to ``best_calc``.

    Returns ``(energies, picks)``: ``energies`` is ``{key: G_eV}`` for the keys
    that have a usable calc (keys without one are left out), and ``picks`` is
    ``{key: Best}`` for every key, so the chosen calc or the reason there is
    none ("no finished calcs", ...) stays inspectable.
    """
    picks = {
        key: (best_calc(df, metal, bias, kind="surf", **kwargs) if name is None
              else best_calc(df, metal, bias, name, **kwargs))
        for key, name in species.items()
    }
    energies = {k: b.row["G_eV"] for k, b in picks.items() if b.row is not None}
    return energies, picks


def candidates(df: pd.DataFrame, metal: str, bias: str, species: str | None = None,
               kind: str = "ads", states: tuple[str, ...] = USABLE_STATES,
               by: str = "counts_as", exclusive: bool = False) -> pd.DataFrame:
    """Every calc :func:`best_calc` considers for (metal, bias, species), with
    why each did or did not qualify. Same arguments as ``best_calc``.

    Columns: ``rank`` (1 = what ``best_calc`` returns; NaN if it failed the
    gate), ``config``, ``state``, ``counts_as``, ``G_eV``, ``dG_vs_best``
    (eV above the winner), ``why_not`` (None if it passed), ``rel``. Passing
    calcs come first, lowest ``G_eV`` first (ties by ``rel``, as in
    ``best_calc``), then the failures by ``G_eV``.
    """
    base = (df["kind"] == kind) & (df["metal"] == metal) & (df["bias"] == bias)
    if kind == "ads":
        if species is None:
            raise ValueError("species is required for kind='ads'")
        base &= (df["adsorbate"] == species) | species_match(df, species, by)
        if exclusive:
            base &= species_match(df, species, by, exclusive)
    out = df[base].copy()
    out["why_not"] = gate_reasons(out, states)
    out = out.sort_values(["G_eV", "rel"])
    out = pd.concat([out[out["why_not"].isna()], out[out["why_not"].notna()]])
    ok = out["why_not"].isna()
    out["rank"] = pd.Series(range(1, ok.sum() + 1), index=out.index[ok], dtype=float)
    out["dG_vs_best"] = out["G_eV"] - out.loc[ok, "G_eV"].min() if ok.any() else float("nan")
    return out[["rank", "config", "state", "counts_as", "G_eV", "dG_vs_best",
                "why_not", "rel"]].reset_index(drop=True)
