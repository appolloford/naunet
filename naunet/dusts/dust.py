from __future__ import annotations
from typing import TYPE_CHECKING
from ..species import Species
from ..reactions.reactiontype import ReactionType

if TYPE_CHECKING:
    from ..reactions.reaction import Reaction


class Dust:
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

        gdens = " + ".join(f"y[IDX_{s.alias}]" for s in species)
        if gdens:
            self.locvars = [f"double gdens = {gdens}", *self.locvars]
        else:
            self.varis = {**self.varis, "gdens": None}

        self._sym_av = None
        self._sym_tgas = None
        self._sym_tdust = None
        self._sym_radfield = None
        self._sym_crrate = None
        self._sym_h2form = None

    @property
    def sym_av(self) -> str:
        return self._sym_av

    @sym_av.setter
    def sym_av(self, av: str) -> None:
        if not isinstance(av, str):
            raise TypeError("Visual extinction symbol must be string")

        if not av:
            raise ValueError("Visual extinction symbol must be non-empty string")
        self._sym_av = av

    @property
    def sym_tgas(self) -> str:
        return self._sym_tgas

    @sym_tgas.setter
    def sym_tgas(self, tgas: str) -> None:
        if not isinstance(tgas, str):
            raise TypeError("Gas temperature symbol must be string")

        if not tgas:
            raise ValueError("Gas temperature symbol must be non-empty string")
        self._sym_tgas = tgas

    @property
    def sym_tdust(self) -> str:
        return self._sym_tdust

    @sym_tdust.setter
    def sym_tdust(self, tdust: str) -> None:
        if not isinstance(tdust, str):
            raise TypeError("Dust temperature symbol must be string")

        if not tdust:
            raise ValueError("Dust temperature symbol must be non-empty string")
        self._sym_tdust = tdust

    @property
    def sym_radfield(self) -> str:
        return self._sym_radfield

    @sym_radfield.setter
    def sym_radfield(self, radfield: str) -> None:
        if not isinstance(radfield, str):
            raise TypeError("Radiation field symbol must be string")

        if not radfield:
            raise ValueError("Radiation field symbol must be non-empty string")
        self._sym_radfield = radfield

    @property
    def sym_crrate(self) -> str:
        return self._sym_crrate

    @sym_crrate.setter
    def sym_crrate(self, crrate: str) -> None:
        if not isinstance(crrate, str):
            raise TypeError("Cosmic-ray ionization rate symbol must be string")

        if not crrate:
            raise ValueError(
                "Cosmic-ray ionization rate symbol must be non-empty string"
            )
        self._sym_crrate = crrate

    @property
    def sym_h2form(self) -> str:
        return self._sym_h2form

    @sym_h2form.setter
    def sym_h2form(self, h2form: str) -> None:
        if not isinstance(h2form, str):
            raise TypeError("H2 formation rate symbol must be string")

        if not h2form:
            raise ValueError("H2 formation rate symbol must be non-empty string")
        self._sym_h2form = h2form

    def rate_depletion(self, reac: Reaction) -> str:

        if not self.sym_tgas:
            raise ValueError("Gas temperature symbol is not set properly")

        if reac.reaction_type != ReactionType.GRAIN_FREEZE:
            raise ValueError("The reaction type is not depletion")

        if len(reac.reactants) != 1:
            raise ValueError("Number of reactants in depletion should be 1.")

        spec = reac.reactants[0]
        a = reac.alpha

        rate = " * ".join(
            [
                f"{a} * pi * rG * rG * gdens",
                f"sqrt(8.0 * kerg * {self.sym_tgas}/ (pi*amu*{spec.A}))",
            ]
        )
        return rate

    def rate_thermal_desorption(self, reac: Reaction) -> str:

        if not self.sym_tdust:
            raise ValueError("Dust temperature symbol is not set properly")

        if reac.reaction_type != ReactionType.GRAIN_DESORB_THERMAL:
            raise ValueError("The reaction type is not thermal desorption")

        if len(reac.reactants) != 1:
            raise ValueError("Number of reactants in thermal desoprtion should be 1.")

        return NotImplemented

    def rate_photon_desorption(self, reac: Reaction) -> str:

        if not self.sym_radfield:
            raise ValueError("Radiation field symbol is not set properly")

        if reac.reaction_type != ReactionType.GRAIN_DESORB_PHOTON:
            raise ValueError("The reaction type is not photon desorption")

        if len(reac.reactants) != 1:
            raise ValueError("Number of reactants in photon desoprtion should be 1.")

        return NotImplemented

    def rate_cosmicray_desorption(self, reac: Reaction) -> str:

        if not self.sym_crrate:
            raise ValueError("Cosmic ray ionization rate symbol is not set properly")

        if reac.reaction_type != ReactionType.GRAIN_DESORB_COSMICRAY:
            raise ValueError("The reaction type is not cosmic-ray desorption")

        if len(reac.reactants) != 1:
            raise ValueError("Number of species in cosmic-ray desoprtion should be 1.")

        return NotImplemented

    def rate_electron_capture(self, reac: Reaction) -> str:

        if not self.sym_tgas:
            raise ValueError("Gas temperature symbol is not set properly")

        if reac.reaction_type != ReactionType.GRAIN_ECAPTURE:
            raise ValueError("The reaction type is not electron capture")

        return NotImplemented

    def rate_recombination(self, reac: Reaction) -> str:

        if not self.sym_tgas:
            raise ValueError("Gas temperature symbol is not set properly")

        if reac.reaction_type != ReactionType.GRAIN_RECOMINE:
            raise ValueError("The reaction type is not recombination")

        return NotImplemented

    def rate_surface_twobody(self, reac: Reaction) -> str:

        if not self.sym_tdust:
            raise ValueError("Dust temperature symbol is not set properly")

        if reac.reaction_type != ReactionType.SURFACE_TWOBODY:
            raise ValueError("The reaction type is not surface two-body reaction")

        if len(reac.reactants) != 2:
            raise ValueError(
                "Number of species in two-body surface reaction should be 2."
            )

        return NotImplemented

    def rate_reactive_desorption(self, reac: Reaction) -> str:

        if not self.sym_tdust:
            raise ValueError("Dust temperature symbol is not set properly")

        if reac.reaction_type != ReactionType.GRAIN_DESORB_REACTIVE:
            raise ValueError("The reaction type is not reactive desorption")

        if len(reac.reactants) != 2:
            raise ValueError(
                "Number of species in reactive desorption reaction should be 2."
            )

        return NotImplemented

    def rate_h2_desorption(self, reac: Reaction) -> str:

        if not self.sym_h2form:
            raise ValueError("H2 formation rate symbol is not set properly")

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
