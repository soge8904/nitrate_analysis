"""Free-energy diagrams: given a start state and an end state, find the chain
of intermediate states between them and report each one's energy relative to
the start (0 eV).

This is a from-scratch skeleton. States are plain strings (``"NO3-"``,
``"*NO3"``, ``"*NO3H"``); the reaction network is a graph of which states can
follow which; energies are looked up from a flat dict standing in for a future
real data source. No GC/CHE pricing rules or bias math live here yet -- see
the project plan for what's deferred and why.
"""

from __future__ import annotations

Graph = dict[tuple[str, str], None]  # (from_state, to_state) -> metadata (unused for now)

EXAMPLE_GRAPH: Graph = {
    ("NO3-", "*NO3"): None,
    ("*NO3", "*NO3H"): None,
}

EXAMPLE_ENERGIES: dict[str, float] = {
    "NO3-": -5.00,
    "*NO3": -6.20,
    "*NO3H": -7.10,
}


def _norm(state: str) -> str:
    """Canonical form for state-string comparisons: caps-agnostic."""
    return state.strip().upper()


def find_path(graph: Graph, start: str, end: str) -> list[str]:
    """Return the first simple path from ``start`` to ``end`` as a list of
    states, walking ``graph`` depth-first. State matching is caps-agnostic;
    returned states are the canonical (uppercased) form. Raises ``ValueError``
    if no path exists."""
    start, end = _norm(start), _norm(end)
    graph = {(_norm(frm), _norm(to)): meta for (frm, to), meta in graph.items()}

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


def get_FED(
    metal: str,
    start: str,
    end: str,
    bias: str,
    *,
    graph: Graph = EXAMPLE_GRAPH,
    energies: dict[str, float] = EXAMPLE_ENERGIES,
) -> list[tuple[str, float | None]]:
    """Free energy diagram from ``start`` to ``end``: each state on the path
    paired with its energy relative to ``start`` (0 eV), or ``None`` if that
    state (or the start state itself) has no energy in ``energies`` -- the
    state still appears on the path, it just has no y-value yet.

    ``metal`` and ``bias`` are accepted for the target call shape but not yet
    applied -- see the project plan's deferred decisions.
    """
    path = find_path(graph, start, end)
    energies = {_norm(state): value for state, value in energies.items()}
    e0 = energies.get(path[0])

    result: list[tuple[str, float | None]] = []
    for state in path:
        e = energies.get(state)
        result.append((state, None if e is None or e0 is None else e - e0))
    return result
