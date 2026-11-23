import matplotlib

matplotlib.use("Agg")

import pytest
from functools import partial

from example_data import EXAMPLE_CONSTANTS, EXAMPLE_ENERGIES
from nitrate_analysis import fed
from nitrate_analysis.network import NO3RR_GRAPH
from nitrate_analysis.plotting import plot_FED, shift_route

get_FEDs = partial(fed.get_FEDs, graph=NO3RR_GRAPH, energies=EXAMPLE_ENERGIES,
                   constants=EXAMPLE_CONSTANTS)


def _routes():
    # the fake data has no *O_NO2 energy; add one so both routes are complete
    energies = {**EXAMPLE_ENERGIES, "*O_NO2": -105.0}
    return get_FEDs("Fe", "NO3-", "*NO2", "0V", energies=energies)


def test_shift_initial_is_noop():
    route = [("A", 0.0), ("B", -1.0), ("C", 0.5)]
    assert shift_route(route, "initial") == route


def test_shift_final_puts_last_state_at_zero():
    route = [("A", 0.0), ("B", -1.0), ("C", 0.5)]
    assert shift_route(route, "final") == [("A", -0.5), ("B", -1.5), ("C", 0.0)]


def test_shift_keeps_none_and_rejects_none_reference():
    assert shift_route([("A", 0.0), ("B", None), ("C", 1.0)], "final")[1] == ("B", None)
    with pytest.raises(ValueError, match="no energy"):
        shift_route([("A", 0.0), ("B", None)], "final")


def test_shift_rejects_bad_zero():
    with pytest.raises(ValueError, match="zero"):
        shift_route([("A", 0.0)], "middle")


def test_plot_draws_one_bar_per_state_and_one_line_per_gap_free_step():
    routes = _routes()
    ax = plot_FED(routes)
    n_states = sum(len(r) for r in routes)
    n_steps = sum(len(r) - 1 for r in routes)
    bars = sum(len(c.get_segments()) for c in ax.collections)
    dashed = [ln for ln in ax.lines if ln.get_linestyle() == "--"]
    assert bars == n_states
    assert len(dashed) == n_steps


def test_plot_zero_final_sets_every_route_end_to_zero():
    ax = plot_FED(_routes(), zero="final")
    assert ax.get_ylabel() == "G relative to final state (eV)"
    texts = [t.get_text() for t in ax.texts]
    assert sum(t == "*NO2\n0.00" for t in texts) == len(_routes())


def test_plot_marks_missing_energy_and_breaks_line():
    route = [("NO3-", 0.0), ("*NO3", None), ("*NO3H", 1.0)]
    ax = plot_FED(route)
    assert any("n/a" in t.get_text() for t in ax.texts)
    assert not [ln for ln in ax.lines if ln.get_linestyle() == "--"]


def test_unzeroable_route_warns_and_is_not_drawn():
    routes = get_FEDs("Fe", "NO3-", "*NO2", "0V")  # split route has no energy
    with pytest.warns(UserWarning, match="not drawn"):
        ax = plot_FED(routes, zero="final")
    assert any("n/a" in t.get_text() for t in ax.texts)


def test_plot_uses_supplied_axes():
    import matplotlib.pyplot as plt
    _, ax = plt.subplots()
    assert plot_FED(_routes(), ax=ax) is ax
