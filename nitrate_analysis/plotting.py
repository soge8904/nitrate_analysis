"""Plot free-energy diagrams from :func:`nitrate_analysis.fed.get_FEDs`.

Draws numbers it is handed and knows nothing about calcs, references or the
workbook. The caller decides whether to show or save the figure.
"""

from __future__ import annotations

import warnings

import matplotlib.pyplot as plt

Route = list[tuple[str, "float | None"]]

BAR_HALF_WIDTH = 0.3


def shift_route(route: Route, zero: str = "initial") -> Route:
    """Re-zero one route (a list of ``(state, G)``) at its ``"initial"`` or
    ``"final"`` state. ``"initial"`` is what :func:`fed.get_FEDs` already
    returns, so it is a no-op. Raises ``ValueError`` if the chosen state has no
    energy (nothing sensible to reference against)."""
    if zero not in ("initial", "final"):
        raise ValueError(f"zero must be 'initial' or 'final', not {zero!r}")
    ref = route[0][1] if zero == "initial" else route[-1][1]
    if ref is None:
        raise ValueError(
            f"cannot zero at the {zero} state {route[0 if zero == 'initial' else -1][0]!r}: "
            "it has no energy"
        )
    return [(state, None if g is None else g - ref) for state, g in route]


def plot_FED(
    routes: list[Route],
    *,
    zero: str = "initial",
    route_labels: list[str] | None = None,
    ax=None,
    title: str | None = None,
):
    """Plot one or more routes as a free-energy diagram and return the Axes.

    ``routes`` is the output of ``fed.get_FEDs`` (a single route from
    ``fed.get_FED`` is accepted too). ``zero`` picks the 0 eV reference:
    ``"initial"`` (default, the start of each route) or ``"final"``.

    Each state is a horizontal bar at its energy, placed by its position along
    the route, with dashed connectors. A state with no energy (``None``) leaves
    a gap: the line breaks and ``n/a`` is written at its position. A route that
    cannot be zeroed (its reference state has no energy) is not drawn: a warning
    is issued and all of its states show ``n/a``. ``ax`` lets
    the caller supply an Axes (for multi-panel figures); otherwise a new one is
    made.
    """
    if routes and isinstance(routes[0], tuple):  # a single route
        routes = [routes]
    shifted = []
    for r in routes:
        try:
            shifted.append(shift_route(r, zero))
        except ValueError as err:  # e.g. zero="final" but the end state has no energy
            warnings.warn(f"route not drawn: {err}")
            shifted.append([(state, None) for state, _ in r])
    routes = shifted
    if route_labels is None:
        route_labels = [" → ".join(s for s, _ in r) for r in routes]
    if len(route_labels) != len(routes):
        raise ValueError("route_labels must have one entry per route")

    if ax is None:
        _, ax = plt.subplots(figsize=(1.8 * max(len(r) for r in routes) + 1.5, 4.5))

    gaps = []
    for i, (route, label) in enumerate(zip(routes, route_labels)):
        color = f"C{i}"
        ax.plot([], [], color=color, label=label)  # legend entry

        last = None  # (x, g) of the previous state if it had an energy
        for x, (state, g) in enumerate(route):
            if g is None:
                gaps.append((x, state, color))
                last = None
                continue
            ax.hlines(g, x - BAR_HALF_WIDTH, x + BAR_HALF_WIDTH, color=color, lw=2.5)
            ax.annotate(f"{state}\n{g:.2f}", (x, g), xytext=(0, 5),
                        textcoords="offset points", ha="center", va="bottom",
                        fontsize=8, color=color)
            if last is not None:
                ax.plot([last[0] + BAR_HALF_WIDTH, x - BAR_HALF_WIDTH], [last[1], g],
                        color=color, ls="--", lw=1)
            last = (x, g)

    ax.margins(x=0.1, y=0.2)
    for x, state, color in gaps:  # after margins so the y-limits are final
        ax.annotate(f"{state}\nn/a", (x, ax.get_ylim()[0]), xytext=(0, 5),
                    textcoords="offset points", ha="center", va="bottom",
                    fontsize=8, color=color)

    ax.axhline(0, color="0.7", lw=0.8, zorder=0)
    ax.set_xticks(range(max(len(r) for r in routes)))
    ax.set_xticklabels([])
    ax.set_xlabel("reaction coordinate")
    ax.set_ylabel(f"G relative to {zero} state (eV)")
    if title:
        ax.set_title(title)
    if len(routes) > 1:
        ax.legend(fontsize=8, frameon=False)
    return ax
