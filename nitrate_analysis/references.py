"""Module for loading and storing fixed thermodynamic reference quantities."""

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd


@dataclass(frozen=True) #frozen keeps from reassigning values to the refs.
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
        #error for no file found:
        if not path.exists():
            raise FileNotFoundError(f"Reference workbook not found: {path}")
        
        try:
            table = pd.read_excel(path, sheet_name=sheet_name)
        except ValueError as exc:
            raise ReferenceError(
                f"Could not read sheet {sheet_name!r} from {path}."
            ) from exc

        required = {"name", "value_eV"}
        missing = required - set(table.columns)
        if missing:
            raise ReferenceError(
                f"{sheet_name!r} is missing required columns: "
                + ", ".join(sorted(missing))
            )

        table = table.dropna(how="all").copy()
        table["name"] = table["name"].astype(str).str.strip() #makes sure all numbers in mol names are strings        

        #check for duplicate refs
        if table["name"].duplicated().any():
            duplicates = sorted(table.loc[table["name"].duplicated(), "name"].unique())
            raise ReferenceError(f"Duplicate reference names: {duplicates}")

        values: dict[str, float] = {}
        metadata: dict[str, dict[str, Any]] = {}


        for _, row in table.iterrows():
            name = row["name"]
            value = row["value_eV"]
            #if value is empty or not a number, raise error
            ##Check if blank, float(value) would not raise an error for NaN.
            if pd.isna(value):
                raise ReferenceError(f"Reference {name!r} has no numeric value.")
            try:
                values[name] = float(value)
            except (ValueError, TypeError):
                raise ReferenceError(f"Reference {name!r} has non-numeric value {value!r}.")

            metadata[name] = {
                column: row[column]
                for column in table.columns
                if column not in {"name", "value_eV"} and not pd.isna(row[column])
            }

        return cls(values=values, metadata=metadata)

    def get(self, name: str) -> float:
        """Return a reference value by name with a useful error message."""
        try:
            return self.values[name]
        except KeyError as exc:
            available = ", ".join(sorted(self.values))
            raise ReferenceError(
                f"Reference {name!r} is not defined. Available references: {available}"
            ) from exc
        
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
    def dG_rxn_redox(self) -> float:
        return self.get("dG_rxn_redox")

    @property
    def dG_rxn_pka_henry(self) -> float:
        return self.get("dG_rxn_pka_henry")

    @property
    def dG_rxn_henry_diss(self) -> float:
        return self.get("dG_rxn_henry_diss")

    def get_nitrate_association_dGrxn(self, route: str) -> dict[str, float]:
        """Return dG_rxn for the reaction: NO3- (aq) + H+(aq) -> HNO3 (g)"""
        
        route = route.lower().strip()

        if route == "redox":
            dG_rxn = self.dG_rxn_redox
        elif route == "pka_henry":
            dG_rxn = self.dG_rxn_pka_henry
        elif route == "henry_diss":
            dG_rxn = self.dG_rxn_henry_diss
        else:
            raise ReferenceError(
                f"Unknown nitrate route {route!r}. Use 'redox' or 'pka_henry' or 'henry_diss'."
            )

        return dG_rxn
