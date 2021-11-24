from __future__ import annotations
from .species import Species


class ThermalProcess:
    def __init__(
        self,
        reactants: list[str],
        rate: str,
        consts: dict = None,
        globs: dict = None,
        varis: dict = None,
        user_var: list = None,
    ) -> None:

        self._reactants = reactants.copy()
        self._rate = rate
        self.consts = consts.copy() if consts else {}
        self.globs = globs.copy() if globs else {}
        self.varis = varis.copy() if varis else {}
        self.user_var = user_var.copy() if user_var else []

    @property
    def reactant_names(self):
        return self._reactants

    @property
    def reactants(self):
        return [Species(r) for r in self._reactants]

    def rate_func(self):
        return self._rate


HCollisionalIonizationCooling = ThermalProcess(
    ["H", "e-"],
    "1.27e-21 * sqrt(y[IDX_TGAS]) / (1.0 + sqrt(y[IDX_TGAS]/1e5)) * exp(-1.578091e5/y[IDX_TGAS])",
)


supported_heating_process = {}

supported_cooling_process = {
    "CIC_H": HCollisionalIonizationCooling,
}


def get_allowed_heating(species: list[Species]) -> dict[str, ThermalProcess]:
    """
    Get a dict of allowed heating processes of input species list from
    all support heating models.

    Args:
        species (list[Species]): A list of species

    Returns:
        dict[str, ThermalProcess]: allowed heating processes.
                                   Dictionary of {name: <heating process>}
    """

    allowed_heating = {}

    return allowed_heating


def get_allowed_cooling(species: list[Species]) -> dict[str, ThermalProcess]:
    """
    Get a dict of allowed cooling processes of input species list from
    all support cooling models.

    Args:
        species (list[Species]): A list of species

    Returns:
        dict[str, ThermalProcess]: allowed cooling processes.
                                   Dictionary of {name: <cooling process>}
    """

    specnamelist = [s.name for s in species]
    allowed_cooling = {}

    for key, value in supported_cooling_process.items():

        if all(r in specnamelist for r in value.reactant_names):

            allowed_cooling[key] = value

    return allowed_cooling
