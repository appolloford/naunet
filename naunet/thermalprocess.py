from __future__ import annotations
from .species import Species

supported_heating_process = {}

supported_cooling_process = {}


class HeatingProcess:
    def __init__(
        self,
        reactants: list[str],
        rate: str,
        consts: dict = None,
        globs: dict = None,
        varis: dict = None,
        user_var: list = None,
    ) -> None:

        self.reactants = [Species(r) for r in reactants]
        self._rate = rate
        self.consts = consts.copy() if consts else {}
        self.globs = globs.copy() if globs else {}
        self.varis = varis.copy() if varis else {}
        self.user_var = user_var.copy() if user_var else []

    def rate_func(self):
        return self._rate


class CoolingProcess:
    def __init__(
        self,
        reactants: list[str],
        rate: str,
        consts: dict = None,
        globs: dict = None,
        varis: dict = None,
        user_var: list = None,
    ) -> None:

        self.reactants = [Species(r) for r in reactants]
        self._rate = rate
        self.consts = consts.copy() if consts else {}
        self.globs = globs.copy() if globs else {}
        self.varis = varis.copy() if varis else {}
        self.user_var = user_var.copy() if user_var else []

    def rate_func(self):
        return self._rate


known_elements = Species.known_elements()

if "H" in known_elements and "e" in known_elements:
    HCollisionalIonizationCooling = CoolingProcess(
        ["H", "e-"],
        "1.27e-21 * sqrt(y[IDX_TGAS]) / (1.0 + sqrt(y[IDX_TGAS]/1e5)) * exp(-1.578091e5/y[IDX_TGAS])",
    )
    supported_cooling_process["CIC_H"] = HCollisionalIonizationCooling
