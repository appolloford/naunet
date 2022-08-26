from __future__ import annotations
from .component import Component, VariableType as vt
from .species import Species


class ThermalProcess(Component):
    def __init__(
        self,
        reactants: list[str | Species],
        rate: str,
    ) -> None:

        super().__init__()

        self._reactants = [
            self._create_species(r) for r in reactants if self._create_species(r)
        ]

        self.temp_min = -1.0
        self.temp_max = -1.0
        self._rate = rate

        self.register("mean_molecular_weight", ("mu", -1.0, vt.param))
        self.register("heat_capacity_ratio", ("gamma", -1.0, vt.param))

    @property
    def reactants(self):
        return self._reactants

    def rateexpr(self):
        return self._rate


HICollisionalIonizationCooling = ThermalProcess(
    ["H", "e-"],
    "1.27e-21 * sqrt(y[IDX_TGAS]) / (1.0 + sqrt(y[IDX_TGAS]/1e5)) * exp(-1.578091e5/y[IDX_TGAS])",
)

HeICollisionalIonizationCooling = ThermalProcess(
    ["He", "e-"],
    "9.38e-22 * sqrt(y[IDX_TGAS]) / (1.0 + sqrt(y[IDX_TGAS]/1e5)) * exp(-2.853354e5/y[IDX_TGAS])",
)

HeIICollisionalIonizationCooling = ThermalProcess(
    ["He+", "e-"],
    "4.95e-22 * sqrt(y[IDX_TGAS]) / (1.0 + sqrt(y[IDX_TGAS]/1e5)) * exp(-6.31515e5/y[IDX_TGAS])",
)

He_2SCollisionalIonizationCooling = ThermalProcess(
    ["He+", "e-", "e-"],
    "5.01e-27 * pow(y[IDX_TGAS], -0.1687) / (1.0 + sqrt(y[IDX_TGAS]/1e5)) * exp(-5.5338e4/y[IDX_TGAS])",
)

HIIRecombinationCooling = ThermalProcess(
    ["H+", "e-"],
    "8.7e-27 * sqrt(y[IDX_TGAS]) * pow(y[IDX_TGAS]/1e3, -0.2) / (1.0+pow(y[IDX_TGAS]/1e6, 0.7))",
)

# Dielectronic recombination cooling
HeIRecombinationCooling = ThermalProcess(
    ["He+", "e-"],
    "1.24e-13 * pow(y[IDX_TGAS], -1.5) * exp(-4.7e5/y[IDX_TGAS]) * (1.0+0.3*exp(-9.4e4/y[IDX_TGAS]))",
)

HeIIRecombinationCooling = ThermalProcess(
    ["He+", "e-"], "1.55e-26 * pow(y[IDX_TGAS], 0.3647)"
)

HeIIIRecombinationCooling = ThermalProcess(
    ["He++", "e-"],
    "3.48e-26 * sqrt(y[IDX_TGAS]) * pow(y[IDX_TGAS]/1e3, -0.2) / (1.0+pow(y[IDX_TGAS]/1e6, 0.7))",
)

HICollisionalExcitationCooling = ThermalProcess(
    ["H", "e-", "e-"],
    "9.1e-27 * pow(y[IDX_TGAS], -0.1687) / (1.0+sqrt(y[IDX_TGAS]/1e5)) * exp(-1.3179e4/y[IDX_TGAS])",
)

HeICollisionalExcitationCooling = ThermalProcess(
    ["He+", "e-"],
    "5.54e-17 * pow(y[IDX_TGAS], -.0397) / (1.0+sqrt(y[IDX_TGAS]/1e5)) *exp(-4.73638e5/y[IDX_TGAS])",
)

HeIICollisionalExcitationCooling = ThermalProcess(
    ["He+", "e-"],
    "5.54e-17 * pow(y[IDX_TGAS], -.0397) / (1.0+sqrt(y[IDX_TGAS]/1e5)) *exp(-4.73638e5/y[IDX_TGAS])",
)

supported_heating_process = {}

supported_cooling_process = {
    "CIC_HI": HICollisionalIonizationCooling,
    "CIC_HeI": HeICollisionalIonizationCooling,
    "CIC_HeII": HeIICollisionalIonizationCooling,
    "CIC_He_2S": He_2SCollisionalIonizationCooling,
    "RC_HII": HIIRecombinationCooling,
    "RC_HeI": HeIRecombinationCooling,
    "RC_HeII": HeIIRecombinationCooling,
    "RC_HeIII": HeIIIRecombinationCooling,
    "CEC_HI": HICollisionalExcitationCooling,
    "CEC_HeI": HeICollisionalExcitationCooling,
    "CEC_HeII": HeIICollisionalExcitationCooling,
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

    allowed_cooling = {}

    for key, value in supported_cooling_process.items():

        if all(r in species for r in value.reactants):

            allowed_cooling[key] = value

    return allowed_cooling
