import pytest

from nitrate_analysis.fed import EXAMPLE_ENERGIES, EXAMPLE_GRAPH, find_path, get_FED


def test_find_path_returns_the_chain_through_intermediates():
    assert find_path(EXAMPLE_GRAPH, "NO3-", "*NO3H") == ["NO3-", "*NO3", "*NO3H"]


def test_find_path_start_equals_end_is_trivial():
    assert find_path(EXAMPLE_GRAPH, "*NO3", "*NO3") == ["*NO3"]


def test_find_path_raises_when_disconnected():
    with pytest.raises(ValueError, match="no path"):
        find_path(EXAMPLE_GRAPH, "*NO3H", "NO3-")
        #need to add catalytic cycle understanding


def test_find_path_follows_a_fork():
    graph = {
        ("A", "B"): None,
        ("A", "C"): None,
        ("C", "D"): None,
    }
    path = find_path(graph, "A", "D")
    assert path == ["A", "C", "D"]


def test_get_fed_references_energies_to_start():
    result = get_FED("Fe", "NO3-", "*NO3H", "0V")
    assert result == [
        ("NO3-", 0.0),
        ("*NO3", pytest.approx(EXAMPLE_ENERGIES["*NO3"] - EXAMPLE_ENERGIES["NO3-"])),
        ("*NO3H", pytest.approx(EXAMPLE_ENERGIES["*NO3H"] - EXAMPLE_ENERGIES["NO3-"])),
    ]
    assert result[0][1] == 0.0


def test_get_fed_uses_injected_graph_and_energies():
    graph = {("X", "Y"): None}
    energies = {"X": -1.0, "Y": -3.5}
    result = get_FED("Fe", "X", "Y", "0V", graph=graph, energies=energies)
    assert result == [("X", 0.0), ("Y", -2.5)]

def test_get_fed_with_intermediate_missing():
    graph = {
        ("A", "B"): None,
        ("B", "C"): None,
    }
    energies = {"A": -5.00, "C": -7.00}
    result = get_FED("Fe", "A", "C", "0V", graph=graph, energies=energies)
    assert result == [("A", 0.0), ("B", None), ("C", pytest.approx(-2.0))]


def test_get_fed_with_start_missing_leaves_everything_none():
    graph = {("A", "B"): None}
    energies = {"B": -3.0}
    result = get_FED("Fe", "A", "B", "0V", graph=graph, energies=energies)
    assert result == [("A", None), ("B", None)]


def test_find_path_is_caps_agnostic():
    assert find_path(EXAMPLE_GRAPH, "no3-", "*no3h") == ["NO3-", "*NO3", "*NO3H"]


def test_get_fed_is_caps_agnostic():
    result = get_FED("Fe", "no3-", "*no3", "0V")
    assert result == [("NO3-", 0.0), ("*NO3", pytest.approx(-1.2))]
