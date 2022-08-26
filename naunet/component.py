from __future__ import annotations
import logging
from collections import OrderedDict
from dataclasses import dataclass
from enum import Enum
from types import SimpleNamespace
from .species import Species


class VariableType(Enum):
    constant = 1
    param = 2
    derived = 3


@dataclass
class Variable:
    symbol: str
    value: str | float | None
    type: VariableType


class Component:
    def __init__(self) -> None:
        self._symbols = OrderedDict()

    def _create_species(self, species_name: Species | str) -> Species:
        """
        Create a Species instance if the name is not a pseudo element
        (e.g. CR, CRPHOT), else return None

        Args:
            species_name (Species | str): name of species

        Returns:
            Species: Species instance of the input name
        """

        if isinstance(species_name, Species):
            return species_name

        if species_name and species_name not in Species.known_pseudoelements():
            return Species(species_name)

    def register(
        self,
        name: str,
        variable: tuple[str, str | float | None, VariableType],
        force_overwrite=False,
    ) -> None:

        if self._symbols.get(name) is None or force_overwrite:
            self._symbols[name] = Variable(*variable)

        else:
            logging.warning(
                f"`{name}` has been registered. "
                "Set `force_overwritten = True` to overwrite."
            )

    @property
    def constants(self) -> dict[str, float]:
        symbols = [
            sym for _, sym in self._symbols.items() if sym.type == VariableType.constant
        ]
        return {sym.symbol: sym.value for sym in symbols}

    @property
    def symbols(self) -> SimpleNamespace:
        return SimpleNamespace(**self._symbols)

    @property
    def params(self) -> dict[str, float]:
        symbols = [
            sym for _, sym in self._symbols.items() if sym.type == VariableType.param
        ]
        return {sym.symbol: sym.value for sym in symbols}

    @property
    def deriveds(self) -> OrderedDict[str, str]:
        symbols = [
            sym for _, sym in self._symbols.items() if sym.type == VariableType.derived
        ]
        return OrderedDict((sym.symbol, sym.value) for sym in symbols)
