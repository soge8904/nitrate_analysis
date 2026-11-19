"""Valence electron counts from a chemical formula.

The counts are set by the pseudopotentials (GBRV), so they are fixed facts, not
results. ``valence_electrons("NO3-")`` = 5 + 3*6 + 1 = 24; ``("H2O")`` = 8.
Formulas are element symbols with optional counts and an optional single
trailing charge sign (``"NO3-"``, ``"NH4+"``). No parentheses or hydrates.
"""

from __future__ import annotations

import re

# Valence electrons per element in the GBRV pseudopotentials. Add elements as
# they are needed (e.g. metals: GBRV Fe includes semicore 3s3p, so 16, not 8).
ZVAL: dict[str, int] = {"H": 1, "N": 5, "O": 6}

_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")


def valence_electrons(formula: str) -> int:
    """Total valence electrons of ``formula``. Raises ``ValueError`` for an
    unknown element or a formula that is not fully understood."""
    body, charge = formula.strip(), 0
    if body.endswith(("+", "-")):
        charge = 1 if body[-1] == "+" else -1
        body = body[:-1]
    tokens = _TOKEN.findall(body)
    if not tokens or "".join(el + n for el, n in tokens) != body:
        raise ValueError(f"cannot parse formula {formula!r}")
    total = 0
    for element, count in tokens:
        if element not in ZVAL:
            raise ValueError(
                f"no valence count for element {element!r} (in {formula!r}); "
                f"known: {', '.join(ZVAL)}"
            )
        total += ZVAL[element] * int(count or 1)
    return total - charge
