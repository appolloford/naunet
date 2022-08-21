from __future__ import annotations
from typing import TYPE_CHECKING
from ..species import Species
from ..reactions.reactiontype import ReactionType

if TYPE_CHECKING:
    from ..reactions.reaction import Reaction


class Grain:
    """
    Class of base dust grain model.

    Attributes:
        model (str): The name of the model
        varis (dict): The parameters exist in the model
        locvars (list): The local variables in the model
    """

    model = "base"

    consts = {}
    varis = {}
    locvars = []

    def __init__(self, species: list[Species] = None, group: int = 0) -> None:

        self.species = species.copy() if species else []

        gdens = " + ".join(f"y[IDX_{s.alias}]" for s in self.species)
        if gdens:
            self.locvars = [f"double gdens = {gdens}", *self.locvars]
        else:
            self.varis = {**self.varis, "gdens": None}

    def rate_depletion(self, reac: Reaction) -> str:

        if reac.reaction_type != ReactionType.GRAIN_FREEZE:
            raise ValueError("The reaction type is not depletion")

        if len(reac.reactants) != 1:
            raise ValueError("Number of reactants in depletion should be 1.")

        spec = reac.reactants[0]
        a = reac.alpha
        tgas = reac.symbols.temperature.symbol

        rate = " * ".join(
            [
                f"{a} * pi * rG * rG * gdens",
                f"sqrt(8.0 * kerg * {tgas}/ (pi*amu*{spec.A}))",
            ]
        )
        return rate

    def rate_thermal_desorption(self, reac: Reaction) -> str:

        if reac.reaction_type != ReactionType.GRAIN_DESORB_THERMAL:
            raise ValueError("The reaction type is not thermal desorption")

        if len(reac.reactants) != 1:
            raise ValueError("Number of reactants in thermal desoprtion should be 1.")

        return NotImplemented

    def rate_photon_desorption(self, reac: Reaction) -> str:

        if reac.reaction_type != ReactionType.GRAIN_DESORB_PHOTON:
            raise ValueError("The reaction type is not photon desorption")

        if len(reac.reactants) != 1:
            raise ValueError("Number of reactants in photon desoprtion should be 1.")

        return NotImplemented

    def rate_cosmicray_desorption(self, reac: Reaction) -> str:

        if reac.reaction_type != ReactionType.GRAIN_DESORB_COSMICRAY:
            raise ValueError("The reaction type is not cosmic-ray desorption")

        if len(reac.reactants) != 1:
            raise ValueError("Number of species in cosmic-ray desoprtion should be 1.")

        return NotImplemented

    def rate_electron_capture(self, reac: Reaction) -> str:

        if reac.reaction_type != ReactionType.GRAIN_ECAPTURE:
            raise ValueError("The reaction type is not electron capture")

        return NotImplemented

    def rate_recombination(self, reac: Reaction) -> str:

        if reac.reaction_type != ReactionType.GRAIN_RECOMINE:
            raise ValueError("The reaction type is not recombination")

        return NotImplemented

    def rate_surface_twobody(self, reac: Reaction) -> str:

        if reac.reaction_type != ReactionType.SURFACE_TWOBODY:
            raise ValueError("The reaction type is not surface two-body reaction")

        if len(reac.reactants) != 2:
            raise ValueError(
                "Number of species in two-body surface reaction should be 2."
            )

        return NotImplemented

    def rate_reactive_desorption(self, reac: Reaction) -> str:

        if reac.reaction_type != ReactionType.GRAIN_DESORB_REACTIVE:
            raise ValueError("The reaction type is not reactive desorption")

        if len(reac.reactants) != 2:
            raise ValueError(
                "Number of species in reactive desorption reaction should be 2."
            )

        return NotImplemented

    def rate_h2_desorption(self, reac: Reaction) -> str:

        if reac.reaction_type != ReactionType.GRAIN_DESORB_H2:
            raise ValueError("The reaction type is not H2 desorption")

        if len(reac.reactants) != 1:
            raise ValueError("Number of species in H2 desoprtion should be 1.")

        return NotImplemented

    def rateexpr(self, reac: Reaction) -> str:

        rtype = reac.reaction_type

        if rtype == ReactionType.GRAIN_RECOMINE:
            rate = self.rate_recombination(reac)

        elif rtype == ReactionType.GRAIN_FREEZE:
            rate = self.rate_depletion(reac)

        elif rtype == ReactionType.GRAIN_DESORB_THERMAL:
            rate = self.rate_thermal_desorption(reac)

        elif rtype == ReactionType.GRAIN_DESORB_PHOTON:
            rate = self.rate_photon_desorption(reac)

        elif rtype == ReactionType.GRAIN_DESORB_COSMICRAY:
            rate = self.rate_cosmicray_desorption(reac)

        elif rtype == ReactionType.GRAIN_DESORB_H2:
            rate = self.rate_h2_desorption(reac)

        elif rtype == ReactionType.SURFACE_TWOBODY:
            rate = self.rate_surface_twobody(reac)

        elif rtype == ReactionType.GRAIN_DESORB_REACTIVE:
            rate = self.rate_reactive_desorption(reac)

        elif rtype == ReactionType.GRAIN_ECAPTURE:
            rate = self.rate_electron_capture(reac)

        else:
            raise ValueError(
                f"Unknown reaction type in {self.model} dust model: {rtype}"
            )

        if rate is NotImplemented:
            raise NotImplementedError(
                f"The reaction rate function is not implemented in {self.model}"
            )

        return rate
