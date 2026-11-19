"""The NO3RR reaction network as data: which states follow which.

States are plain strings: aqueous (``"NO3-"``) or adsorbed (``"*NO3"``, with
``"*"`` the bare slab). Adsorbed names match the calculation species names
(``calcs`` ``counts_as``) after dropping the ``*``. ``fed.py`` walks and prices
a graph; it knows nothing about this particular one.
"""

from nitrate_analysis.fed import Graph

# (from_state, to_state) -> per-step (n_h2o, n_h, and optional "ref"), or
# None for the adsorption-type edge (which needs no per-edge numbers -- see
# nitrate_adsorption_dG). Two routes to *NO2:
#   * intact: NO3- -> *NO3 -> *NO3H -> *NO2. *NO2 is only reachable through
#     *NO3H (no direct *NO3 -> *NO2 edge), but its dG is still priced against
#     *NO3 via "ref".
#   * split: NO3- -> *O_NO2 (adsorbs and splits into *O + *NO2, same atoms as
#     *NO3) -> *NO2, losing the O as water (*O + 2(H+ + e-) -> H2O).
# The first-listed route is the one get_FED / find_path return.
NO3RR_GRAPH: Graph = {
    ("NO3-", "*NO3"): None,
    ("NO3-", "*O_NO2"): None,
    ("*NO3", "*NO3H"): {"n_h2o": 0, "n_h": 1},
    ("*NO3H", "*NO2"): {"n_h2o": 1, "n_h": 2, "ref": "*NO3"},
    ("*O_NO2", "*NO2"): {"n_h2o": 1, "n_h": 2},
}
