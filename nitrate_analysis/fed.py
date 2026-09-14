"""Free-energy diagrams: given a start state and an end state, find the chain
of intermediate states between them and report each one's energy relative to
the start (0 eV).

States are plain strings (``"NO3-"``, ``"*NO3"``, ``"*NO3H"``); the reaction
network is a graph of which states can follow which. Two different equations
price a step, chosen by whether the state names start with ``"*"``
(adsorbed) or not (aqueous):

* aqueous -> adsorbed (currently only ``"NO3-" -> "*NO3"``): the nitrate
  adsorption step, priced by :func:`nitrate_adsorption_dG`.
* adsorbed -> adsorbed (e.g. ``"*NO3" -> "*NO3H"``): a proton/water step,
  priced by :func:`pcet_dG` using a per-edge hardcoded ``(n_h2o, n_h)`` --
  how many waters and protons are transferred relative to the earlier state.
  These are hand-specified per edge (as graph metadata) rather than derived
  from formulas. The electron term is *not* an independent per-edge number:
  it comes from how many valence electrons leave with each transferred
  water (``n_h2o * N_H2O``) -- see ``EXAMPLE_VALENCE_ELECTRONS`` below.

Graph *connectivity* (which state must be visited before which) and a step's
*energy reference point* (which state's G the dG formula subtracts) are two
different things that don't have to agree. ``*NO3 -> *NO3H -> *NO2`` means
you must visit ``*NO3H`` to reach ``*NO2`` -- but ``*NO2``'s dG is still
written against ``*NO3`` directly, not ``*NO3H``. An edge's metadata can
carry an optional ``"ref"`` key naming that state; it defaults to the edge's
own predecessor when omitted (the common case, where they're the same).
:func:`get_FED` looks up the *reference* state's cumulative energy to build
on, not necessarily its immediate predecessor's -- so a state missing its
own energy only blocks whatever actually depends on it, not everything
physically downstream.

Energies (``EXAMPLE_ENERGIES``), universal constants (``EXAMPLE_CONSTANTS``,
e.g. ``mu_H+``, ``mu_H2O``), and tabulated valence-electron counts of
reference molecules (``EXAMPLE_VALENCE_ELECTRONS``, e.g. ``N_NO3- = 24``,
``N_H2O = 8``) are all fake local stand-ins for a future real data source,
kept in separate dicts since they are different kinds of numbers.
"""

from __future__ import annotations

Graph = dict[tuple[str, str], dict[str, float | str] | None]

# (from_state, to_state) -> per-step (n_h2o, n_h, and optional "ref"), or
# None for the adsorption-type edge (which needs no per-edge numbers -- see
# nitrate_adsorption_dG). *NO2 is only reachable through *NO3H (no direct
# *NO3 -> *NO2 edge), but its dG is still priced against *NO3 via "ref".
EXAMPLE_GRAPH: Graph = {
    ("NO3-", "*NO3"): None,
    ("*NO3", "*NO3H"): {"n_h2o": 0, "n_h": 1},
    ("*NO3H", "*NO2"): {"n_h2o": 1, "n_h": 2, "ref": "*NO3"},
}

# Per-state DFT/reference energies (eV). "*" is the bare slab; "NO3-" has no
# entry -- its energy never comes from DFT, only from the constants below.
EXAMPLE_ENERGIES: dict[str, float] = {
    "*": -100.0,
    "*NO3": -106.0,
    "*NO3H": -106.9,
    "*NO2": -104.0,
}

# Universal constants (eV), not tied to a single state.
EXAMPLE_CONSTANTS: dict[str, float] = {
    "mu_e0": -4.5,  # mu_e at U = 0 V; get_FED's bias shifts this
    "mu_H+": -0.5,
    "mu_H2O": -14.0,
    "mu_HNO3": -18.0,
    "dG_HNO3": 0.3,  # NO3-(aq) + H+(aq) -> HNO3(g)
}

# Valence electron counts of reference molecules (fixed facts about the
# molecule, not DFT results) -- looked up wherever an equation needs "N_X".
EXAMPLE_VALENCE_ELECTRONS: dict[str, float] = {
    "NO3-": 24,
    "H2O": 8,
}


def _norm(state: str) -> str:
    """Canonical form for state-string comparisons: caps-agnostic."""
    return state.strip().upper()


def _normalize_graph(graph: Graph) -> Graph:
    return {(_norm(frm), _norm(to)): meta for (frm, to), meta in graph.items()}


def find_path(graph: Graph, start: str, end: str) -> list[str]:
    """Return the first simple path from ``start`` to ``end`` as a list of
    states, walking ``graph`` depth-first. State matching is caps-agnostic;
    returned states are the canonical (uppercased) form. Raises ``ValueError``
    if no path exists."""
    start, end = _norm(start), _norm(end)
    graph = _normalize_graph(graph)

    if start == end:
        return [start]

    def dfs(node: str, visited: set[str]) -> list[str] | None:
        for frm, to in graph:
            if frm != node or to in visited:
                continue
            if to == end:
                return [node, to]
            found = dfs(to, visited | {to})
            if found is not None:
                return [node, *found]
        return None

    path = dfs(start, {start})
    if path is None:
        raise ValueError(f"no path from {start!r} to {end!r} in graph")
    return path


def _parse_bias(bias: str) -> float:
    return float(bias.strip().rstrip("Vv"))


def _mu_e(constants: dict[str, float], bias: str) -> float | None:
    mu_e0 = constants.get("mu_e0")
    if mu_e0 is None:
        return None
    return mu_e0 - _parse_bias(bias)


def nitrate_adsorption_dG(
    to_state: str,
    bias: str,
    *,
    energies: dict[str, float] = EXAMPLE_ENERGIES,
    constants: dict[str, float] = EXAMPLE_CONSTANTS,
    valence_electrons: dict[str, float] = EXAMPLE_VALENCE_ELECTRONS,
) -> float | None:
    """``NO3-(aq) + * -> to_state + e-``:
    ``dG = G_to - G_* - mu_e*(-N_NO3-) - dG_HNO3 + mu_H+ - mu_HNO3``, where
    ``N_NO3-`` is nitrate's tabulated valence electron count (24). Specific
    to nitrate adsorption -- not a general aqueous-adsorption formula.
    Returns ``None`` if any required value is missing, rather than raising.
    """
    g_bare = energies.get("*")
    g_to = energies.get(_norm(to_state))
    mu_e = _mu_e(constants, bias)
    mu_h = constants.get("mu_H+")
    mu_hno3 = constants.get("mu_HNO3")
    dg_hno3 = constants.get("dG_HNO3")
    n_no3 = valence_electrons.get("NO3-")
    if None in (g_bare, g_to, mu_e, mu_h, mu_hno3, dg_hno3, n_no3):
        return None
    return g_to - g_bare - mu_e * (-n_no3) - dg_hno3 + mu_h - mu_hno3


def pcet_dG(
    frm: str,
    to: str,
    metadata: dict[str, float | str] | None,
    bias: str,
    *,
    energies: dict[str, float] = EXAMPLE_ENERGIES,
    constants: dict[str, float] = EXAMPLE_CONSTANTS,
    valence_electrons: dict[str, float] = EXAMPLE_VALENCE_ELECTRONS,
) -> float | None:
    """``*frm + n_h (H+) -> *to + n_h2o H2O`` (n_h2o/n_h taken from
    ``metadata``, hardcoded per edge):
    ``dG = G_to - G_frm - n_h2o*N_H2O*mu_e - n_h*mu_H+ + n_h2o*mu_H2O``.

    The electron term is driven entirely by the water transferred -- there is
    no independent electron count. A step with no water (``n_h2o=0``, e.g. a
    bare protonation) has no electron term at all.

    Returns ``None`` if ``metadata`` is ``None`` or any required value is
    missing, rather than raising.
    """
    if metadata is None:
        return None
    g_frm = energies.get(_norm(frm))
    g_to = energies.get(_norm(to))
    mu_e = _mu_e(constants, bias)
    mu_h = constants.get("mu_H+")
    mu_h2o = constants.get("mu_H2O")
    n_h2o_electrons = valence_electrons.get("H2O")
    if None in (g_frm, g_to, mu_e, mu_h, mu_h2o, n_h2o_electrons):
        return None
    n_h2o = metadata["n_h2o"]
    return (
        g_to - g_frm
        - n_h2o * n_h2o_electrons * mu_e
        - metadata["n_h"] * mu_h
        + n_h2o * mu_h2o
    )


def get_FED(
    metal: str,
    start: str,
    end: str,
    bias: str,
    *,
    graph: Graph = EXAMPLE_GRAPH,
    energies: dict[str, float] = EXAMPLE_ENERGIES,
    constants: dict[str, float] = EXAMPLE_CONSTANTS,
    valence_electrons: dict[str, float] = EXAMPLE_VALENCE_ELECTRONS,
) -> list[tuple[str, float | None]]:
    """Free energy diagram from ``start`` to ``end``: each state on the path
    paired with its energy relative to ``start`` (0 eV).

    Each step's ``dG`` is priced by :func:`nitrate_adsorption_dG` or
    :func:`pcet_dG`, then added to the cumulative energy of that edge's
    *reference* state (its metadata's ``"ref"``, or its own predecessor in
    the path if there isn't one) -- not necessarily the immediately
    preceding state. A state missing its own energy shows as ``None``; that
    only blocks a later state if the later state's edge actually references
    it -- a state that's merely visited on the way, but not depended on,
    doesn't poison anything after it.

    ``metal`` is accepted for the target call shape but not yet applied.
    """
    path = find_path(graph, start, end)
    graph = _normalize_graph(graph)

    # Insertion order matches path order and find_path never revisits a
    # state, so this dict doubles as the ordered (state, energy) result --
    # list(cumulative.items()) at the end needs no separate list alongside it.
    cumulative: dict[str, float | None] = {path[0]: 0.0}
    for frm, to in zip(path, path[1:]):
        metadata = graph.get((frm, to))
        ref = _norm(metadata["ref"]) if metadata and "ref" in metadata else frm
        ref_value = cumulative.get(ref)
        if ref_value is None:
            cumulative[to] = None
            continue
        if to.startswith("*") and not ref.startswith("*"):
            step = nitrate_adsorption_dG(
                to, bias, energies=energies, constants=constants, valence_electrons=valence_electrons
            )
        elif ref.startswith("*") and to.startswith("*"):
            step = pcet_dG(
                ref, to, metadata, bias,
                energies=energies, constants=constants, valence_electrons=valence_electrons,
            )
        else:
            raise ValueError(f"don't know how to price {ref!r} -> {to!r}")
        cumulative[to] = None if step is None else ref_value + step
    return list(cumulative.items())

FED_energies = get_FED("Ag", start="NO3-", end="*NO2", bias="0V")