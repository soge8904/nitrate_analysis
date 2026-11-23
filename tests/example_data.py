"""Fake energies/constants for testing fed.py (not real calculations)."""

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
    "mu_NO3_H": -17.0,  # chemical potential of the (NO3- + H+) pair
    "mu_NO2_H": -12.0,  # chemical potential of the (NO2- + H+) pair
    "dG_rxn_redox": 0.3,  # NO3-(aq) + H+(aq) -> HNO3(g)
}
