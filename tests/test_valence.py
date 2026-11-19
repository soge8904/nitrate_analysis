import pytest

from nitrate_analysis.valence import valence_electrons


@pytest.mark.parametrize("formula, n", [
    ("H2O", 8), ("NO3-", 24), ("NO2-", 18), ("NH4+", 8), ("HNO3", 24), ("N", 5),
])
def test_known_formulas(formula, n):
    assert valence_electrons(formula) == n


def test_unknown_element_raises():
    with pytest.raises(ValueError, match="Xx"):
        valence_electrons("Xx2")


@pytest.mark.parametrize("bad", ["", "no3", "NO3--", "(NO3)2"])
def test_unparseable_raises(bad):
    with pytest.raises(ValueError):
        valence_electrons(bad)
