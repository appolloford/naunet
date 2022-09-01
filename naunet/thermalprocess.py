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
        self.register("thermalprocess_temperature", ("Temp", "y[IDX_TGAS]", vt.derived))
        self.register("temperature_1e3", ("T3", "Temp/1e3", vt.derived))
        self.register("temperature_1e5", ("T5", "Temp/1e5", vt.derived))
        self.register("temperature_1e6", ("T6", "Temp/1e6", vt.derived))

    @property
    def reactants(self):
        return self._reactants

    def rateexpr(self):
        return self._rate


HICollisionalIonizationCooling = ThermalProcess(
    ["H", "e-"],
    "1.27e-21 * sqrt(Temp) / (1.0 + sqrt(T5)) * exp(-1.578091e5/Temp)",
)

HeICollisionalIonizationCooling = ThermalProcess(
    ["He", "e-"],
    "9.38e-22 * sqrt(Temp) / (1.0 + sqrt(T5)) * exp(-2.853354e5/Temp)",
)

HeIICollisionalIonizationCooling = ThermalProcess(
    ["He+", "e-"],
    "4.95e-22 * sqrt(Temp) / (1.0 + sqrt(T5)) * exp(-6.31515e5/Temp)",
)

He_2SCollisionalIonizationCooling = ThermalProcess(
    ["He+", "e-", "e-"],
    "5.01e-27 * pow(Temp, -0.1687) / (1.0 + sqrt(T5)) * exp(-5.5338e4/Temp)",
)

HIIRecombinationCooling = ThermalProcess(
    ["H+", "e-"],
    "8.7e-27 * sqrt(Temp) * pow(T3, -0.2) / (1.0+pow(T6, 0.7))",
)

# Dielectronic recombination cooling
HeIRecombinationCooling = ThermalProcess(
    ["He+", "e-"],
    "1.24e-13 * pow(Temp, -1.5) * exp(-4.7e5/Temp) * (1.0+0.3*exp(-9.4e4/Temp))",
)

HeIIRecombinationCooling = ThermalProcess(["He+", "e-"], "1.55e-26 * pow(Temp, 0.3647)")

HeIIIRecombinationCooling = ThermalProcess(
    ["He++", "e-"],
    "3.48e-26 * sqrt(Temp) * pow(T3, -0.2) / (1.0+pow(T6, 0.7))",
)

HICollisionalExcitationCooling = ThermalProcess(
    ["H", "e-"], "7.5e-19 / (1.0+sqrt(T5)) * exp(-1.18348e5 / Temp)"
)

HeICollisionalExcitationCooling = ThermalProcess(
    ["He+", "e-"],
    "9.1e-27 * pow(Temp, -0.1687) / (1.0+sqrt(T5)) * exp(-1.3179e4/Temp)",
)

HeIICollisionalExcitationCooling = ThermalProcess(
    ["He+", "e-"],
    "5.54e-17 * pow(Temp, -0.397) / (1.0+sqrt(T5)) *exp(-4.73638e5/Temp)",
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
