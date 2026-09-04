"""Module for loading and storing fixed thermodynamic reference quantities."""

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd


@dataclass
class ReferenceData:
    """Container for fixed thermodynamic reference quantities.

    Values are stored in eV per molecule (or eV per relevant reaction unit),
    matching the convention used in the workbook.
    """

    values: dict[str, float]
    metadata: dict[str, dict[str, Any]]

    @classmethod
    def from_excel(
        cls,
        path: str | Path,
        sheet_name: str = "Reference_summary",
    ) -> "ReferenceData":
        """Load reference values from the workbook summary sheet.

        Parameters
        ----------
        path:
            Path to the Excel workbook.
        sheet_name:
            Name of the summary sheet.

        Returns
        -------
        ReferenceData
            Reference values and their metadata.
        """
        path = Path(path)

        values: dict[str, float] = {}
        metadata: dict[str, dict[str, Any]] = {}

        table = table.dropna(how="all").copy()
        table["name"] = table["name"].astype(str).str.strip()

        for _, row in table.iterrows():
            name = row["name"]
            value = row["value_eV"]

            # Blank rows are ignored, but a named row must have a numeric value.
            values[name] = float(value)

            metadata[name] = {
                column: row[column]
                for column in table.columns
                if column not in {"name", "value_eV"} and not pd.isna(row[column])
            }

        return cls(values=values, metadata=metadata)

    def get(self, name: str) -> float:
        """Return a reference value by name with a useful error message."""
        return self.values[name]


    # --- Fixed references used throughout the analysis ---
    @property
    def mu_H2(self) -> float:
        return self.get("mu_H2")

    @property
    def mu_H2O(self) -> float:
        return self.get("mu_H2O")

    @property
    def mu_NH3(self) -> float:
        return self.get("mu_NH3")

    @property
    def mu_NO(self) -> float:
        return self.get("mu_NO")

    @property
    def mu_NO2(self) -> float:
        return self.get("mu_NO2")

    @property
    def mu_HNO3(self) -> float:
        return self.get("mu_HNO3")

    @property
    def dG_rxn_henry_pka(self) -> float:
        return self.get("dG_rxn_henry_pka")

    @property
    def dG_rxn_redox(self) -> float:
        return self.get("dG_rxn_redox")
