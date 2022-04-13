from __future__ import annotations
from ..species import Species
from ..reactions.reactiontype import ReactionType


class Dust:

    consts = {}
    varis = {}
    locvars = []

    def __init__(
        self,
        species: list[str] | list[Species] = None,
        model: str = "",
    ) -> None:

        self.model = model

        # if dust species have been recorded in Species, use it if the species
        # is not provided. Otherwise the input species overwrite the saved list
        species = species if species else Species.dust_species()
        species = [s if isinstance(s, Species) else Species(s) for s in species]
        # sync the dust species list
        Species.set_dust_species([s.name for s in species])

        self.species = species

        gdens = " + ".join(f"y[IDX_{s.alias}]" for s in species)
        if gdens:
            self.locvars = [f"double gdens = {gdens}", *self.locvars]
        else:
            self.varis = {**self.varis, "gdens": None}

    def rateexpr(
        self,
        rtype: ReactionType,
        reactants: list[Species] = None,
        alpha: float = 0.0,
        beta: float = 0.0,
        gamma: float = 0.0,
        sym_tgas: str = "",
        sym_tdust: str = "",
        sym_phot: str = "",
        sym_cr: str = "",
        **kwargs,
    ) -> str:

        return "0.0"
