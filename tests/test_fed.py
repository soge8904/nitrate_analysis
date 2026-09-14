import pytest

from nitrate_analysis.fed import (
    EXAMPLE_CONSTANTS,
    EXAMPLE_ENERGIES,
    EXAMPLE_GRAPH,
    EXAMPLE_VALENCE_ELECTRONS,
    find_path,
    get_FED,
    nitrate_adsorption_dG,
    pcet_dG,
)


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


def test_find_path_is_caps_agnostic():
    assert find_path(EXAMPLE_GRAPH, "no3-", "*no3h") == ["NO3-", "*NO3", "*NO3H"]


def test_find_path_requires_passing_through_no3h_to_reach_no2():
    # There's no direct *NO3 -> *NO2 edge: you must visit *NO3H first.
    assert find_path(EXAMPLE_GRAPH, "NO3-", "*NO2") == ["NO3-", "*NO3", "*NO3H", "*NO2"]


# -- nitrate_adsorption_dG: NO3-(aq) + * -> *NO3 + e- ----------------------

def test_nitrate_adsorption_dG_matches_formula():
    # dG = G_*NO3 - G_* - mu_e*(-N_NO3-) - dG_HNO3 + mu_H+ - mu_HNO3
    n_no3 = EXAMPLE_VALENCE_ELECTRONS["NO3-"]
    expected = (
        EXAMPLE_ENERGIES["*NO3"]
        - EXAMPLE_ENERGIES["*"]
        - EXAMPLE_CONSTANTS["mu_e0"] * (-n_no3)
        - EXAMPLE_CONSTANTS["dG_HNO3"]
        + EXAMPLE_CONSTANTS["mu_H+"]
        - EXAMPLE_CONSTANTS["mu_HNO3"]
    )
    assert nitrate_adsorption_dG("*NO3", "0V") == pytest.approx(expected)


def test_nitrate_adsorption_dG_varies_with_bias():
    # mu_e = mu_e0 - U, and the electron term is -mu_e*(-N_NO3-) = N_NO3-*mu_e,
    # so a 1V bias shift changes dG by -N_NO3- * 1V.
    n_no3 = EXAMPLE_VALENCE_ELECTRONS["NO3-"]
    at_0V = nitrate_adsorption_dG("*NO3", "0V")
    at_1V = nitrate_adsorption_dG("*NO3", "1V")
    assert at_1V == pytest.approx(at_0V - n_no3)


def test_nitrate_adsorption_dG_accepts_bias_without_v_suffix():
    assert nitrate_adsorption_dG("*NO3", "0") == nitrate_adsorption_dG("*NO3", "0V")


def test_nitrate_adsorption_dG_none_when_target_energy_missing():
    assert nitrate_adsorption_dG("*NO3", "0V", energies={"*": -100.0}) is None


def test_nitrate_adsorption_dG_none_when_valence_electrons_missing():
    assert nitrate_adsorption_dG("*NO3", "0V", valence_electrons={}) is None


# -- pcet_dG: *frm + n_h (H+) -> *to + n_h2o H2O ---------------------------

def test_pcet_dG_pure_protonation_has_no_electron_term():
    # *NO3 -> *NO3H adds one H+ and no water, so the electron term
    # (n_h2o * N_H2O * mu_e) vanishes -- only mu_H+ enters.
    metadata = {"n_h2o": 0, "n_h": 1}
    expected = (
        EXAMPLE_ENERGIES["*NO3H"]
        - EXAMPLE_ENERGIES["*NO3"]
        - EXAMPLE_CONSTANTS["mu_H+"]
    )
    assert pcet_dG("*NO3", "*NO3H", metadata, "0V") == pytest.approx(expected)


def test_pcet_dG_water_producing_step():
    # *NO3 -> *NO2 + H2O: 2 H+ transferred, 1 water produced, and the
    # electron term is driven by that water's valence electrons.
    metadata = {"n_h2o": 1, "n_h": 2}
    n_h2o_electrons = EXAMPLE_VALENCE_ELECTRONS["H2O"]
    expected = (
        EXAMPLE_ENERGIES["*NO2"]
        - EXAMPLE_ENERGIES["*NO3"]
        - 1 * n_h2o_electrons * EXAMPLE_CONSTANTS["mu_e0"]
        - 2 * EXAMPLE_CONSTANTS["mu_H+"]
        + EXAMPLE_CONSTANTS["mu_H2O"]
    )
    assert pcet_dG("*NO3", "*NO2", metadata, "0V") == pytest.approx(expected)


def test_pcet_dG_none_without_metadata():
    assert pcet_dG("*NO3", "*NO3H", None, "0V") is None


def test_pcet_dG_none_when_state_energy_missing():
    metadata = {"n_h2o": 0, "n_h": 1}
    assert pcet_dG("*NO3", "*NO3H", metadata, "0V", energies={"*NO3": -106.0}) is None


def test_pcet_dG_none_when_valence_electrons_missing_and_water_involved():
    metadata = {"n_h2o": 1, "n_h": 2}
    assert pcet_dG("*NO3", "*NO2", metadata, "0V", valence_electrons={}) is None


# -- get_FED: full path, referenced back through *NO3 to the true start ---

def test_get_fed_full_path_from_aqueous_start():
    result = get_FED("Fe", "NO3-", "*NO3H", "0V")
    dG_ads = nitrate_adsorption_dG("*NO3", "0V")
    dG_step = pcet_dG("*NO3", "*NO3H", {"n_h2o": 0, "n_h": 1}, "0V")
    assert result == [
        ("NO3-", 0.0),
        ("*NO3", pytest.approx(dG_ads)),
        ("*NO3H", pytest.approx(dG_ads + dG_step)),
    ]


def test_get_fed_is_caps_agnostic():
    result = get_FED("Fe", "no3-", "*no3", "0V")
    assert result == [("NO3-", 0.0), ("*NO3", pytest.approx(nitrate_adsorption_dG("*NO3", "0V")))]


def test_get_fed_uses_injected_graph_energies_and_constants():
    graph = {("*X", "*Y"): {"n_h2o": 0, "n_h": 1}}
    energies = {"*X": -10.0, "*Y": -11.5}
    constants = {"mu_e0": -4.0, "mu_H+": -0.5, "mu_H2O": -14.0, "mu_HNO3": -18.0, "dG_HNO3": 0.3}
    result = get_FED(
        "Fe", "*X", "*Y", "0V", graph=graph, energies=energies, constants=constants
    )
    # dG = -11.5 - (-10.0) - 0 (no water, no electron term) - 1*(-0.5) = -1.0
    assert result == [("*X", 0.0), ("*Y", pytest.approx(-1.0))]


def test_get_fed_propagates_none_when_energy_missing():
    energies = {"*": -100.0, "*NO3": -106.0}  # *NO3H missing
    result = get_FED("Fe", "*NO3", "*NO3H", "0V", energies=energies)
    assert result == [("*NO3", 0.0), ("*NO3H", None)]


def test_get_fed_propagates_none_when_metadata_missing():
    graph = {("*NO3", "*NO3H"): None}
    result = get_FED("Fe", "*NO3", "*NO3H", "0V", graph=graph)
    assert result == [("*NO3", 0.0), ("*NO3H", None)]


def test_get_fed_propagates_none_past_a_broken_step():
    # A None in the middle of the path must stay None for every later state
    # that actually depends on it (here, *NO2's edge has no "ref", so it
    # defaults to its literal predecessor *NO3H).
    graph = {("*NO3", "*NO3H"): None, ("*NO3H", "*NO2"): {"n_h2o": 0, "n_h": -1}}
    result = get_FED("Fe", "*NO3", "*NO2", "0V", graph=graph)
    assert result == [("*NO3", 0.0), ("*NO3H", None), ("*NO2", None)]


# -- "ref": graph connectivity vs. a step's energy reference state ---------

def test_get_fed_no2_is_priced_against_no3_not_no3h():
    # *NO2's edge carries ref="*NO3" in EXAMPLE_GRAPH, so its dG is computed
    # against *NO3 directly, exactly matching a standalone pcet_dG("*NO3",
    # "*NO2", ...) call -- it is not built by chaining through *NO3H's dG.
    result = get_FED("Fe", "NO3-", "*NO2", "0V")
    dG_ads = nitrate_adsorption_dG("*NO3", "0V")
    dG_no2 = pcet_dG("*NO3", "*NO2", {"n_h2o": 1, "n_h": 2}, "0V")
    assert result[-1] == ("*NO2", pytest.approx(dG_ads + dG_no2))


def test_get_fed_ref_lets_a_later_state_skip_a_missing_predecessor():
    # *NO3H sits on the path (you must visit it to reach *NO2), but its own
    # energy is missing here. Because *NO2's edge references *NO3 (not
    # *NO3H), *NO2 still resolves to a real number -- only *NO3H shows None.
    energies = dict(EXAMPLE_ENERGIES)
    del energies["*NO3H"]
    result = get_FED("Fe", "NO3-", "*NO2", "0V", energies=energies)
    dG_ads = nitrate_adsorption_dG("*NO3", "0V")
    dG_no2 = pcet_dG("*NO3", "*NO2", {"n_h2o": 1, "n_h": 2}, "0V")
    assert result == [
        ("NO3-", 0.0),
        ("*NO3", pytest.approx(dG_ads)),
        ("*NO3H", None),
        ("*NO2", pytest.approx(dG_ads + dG_no2)),
    ]


def test_get_fed_ref_defaults_to_the_literal_predecessor_when_absent():
    # Minimal, graph-agnostic demonstration of the "ref" mechanism itself:
    # *B has no energy, but *C's edge names "ref": "*A", so *C is priced
    # against *A directly and isn't blocked by *B's gap. (Names need the "*"
    # prefix so _step_dG dispatches them as adsorbed -> adsorbed steps.)
    graph = {("*A", "*B"): {"n_h2o": 0, "n_h": 1}, ("*B", "*C"): {"n_h2o": 0, "n_h": 1, "ref": "*A"}}
    energies = {"*A": -10.0, "*C": -8.5}  # *B deliberately has no energy
    constants = {"mu_e0": -4.0, "mu_H+": -0.5, "mu_H2O": -14.0, "mu_HNO3": -18.0, "dG_HNO3": 0.3}
    result = get_FED("Fe", "*A", "*C", "0V", graph=graph, energies=energies, constants=constants)
    # dG(*A->*C) = G_*C - G_*A - 1*mu_H+ = -8.5 - (-10.0) - (-0.5) = 2.0
    assert result == [("*A", 0.0), ("*B", None), ("*C", pytest.approx(2.0))]
