from __future__ import annotations
from .dust import Dust
from ..species import Species
from ..reactions.reactiontype import ReactionType


class RR07Dust(Dust):

    varis = {
        "rG": 1e-5,  # grain radius
        "gdens": 7.6394373e-13,  # grain density
        "sites": 1e15,  # surface sites
        "fr": 1.0,  # freeze ratio
        "opt_thd": 1.0,  # thermal desorption option
        "opt_crd": 1.0,  # cosmic ray induced desorption option
        "opt_h2d": 1.0,  # H2 formation induced desorption option
        "opt_uvd": 1.0,  # UV desorption option
        "eb_h2d": 1.21e3,  # maxmium binding energy which H2 desorption can desorb
        "eb_crd": 1.21e3,  # maxmium binding energy which CR desorption can desorb
        "eb_uvd": 1.0e4,  # maxmium binding energy which UV desorption can desorb
        "crdeseff": 1e5,  # cosmic ray desorption efficiency
        "h2deseff": 1.0e-2,  # H2 desorption efficiency
    }
    locvars = [
        # "double mant = GetMantleDens(y) > 0.0 ? GetMantleDens(y) : 1e-40",
        "double mant = GetMantleDens(y)",
        "double mantabund = mant / nH",
        "double garea = (pi*rG*rG) * gdens",  # total grain cross-section
        "double garea_per_H = garea / nH",
        "double densites = 4.0 * garea * sites",
    ]

    def __init__(self, species: list[str] | list[Species] = None) -> None:
        super().__init__(model="rr07", species=species)

    def _rate_depletion(
        self,
        spec: list[Species],
        a: float,
        b: float,
        c: float,
        sym_tgas: str,
    ) -> str:

        if not isinstance(sym_tgas, str):
            raise TypeError("sym_tgas should be a string")

        if not sym_tgas:
            raise ValueError("Symbol of gas temperature does not exist")

        if len(spec) != 1:
            raise ValueError("Number of species in depletion should be 1.")

        spec = spec[0]
        if spec.iselectron:
            rate = " * ".join(
                [
                    f"4.57e4 * {a} * garea * fr",
                    f"( 1.0 + 16.71e-4/(rG * {sym_tgas}) )",
                ]
            )
        elif spec.charge == 0:
            rate = " * ".join(
                [
                    f"4.57e4 * {a} * garea * fr",
                    f"sqrt({sym_tgas} / {spec.massnumber})",
                ]
            )
        else:
            rate = " * ".join(
                [
                    f"4.57e4 * {a} * garea * fr",
                    f"sqrt({sym_tgas} / {spec.massnumber})",
                    f"( 1.0 + 16.71e-4/(rG * {sym_tgas}) )",
                ]
            )

        return rate

    def _rate_thermal_desorption(
        self,
        spec: list[Species],
        a: float,
        b: float,
        c: float,
        sym_tdust: str,
    ) -> str:

        if not isinstance(sym_tdust, str):
            raise TypeError("sym_tdust should be a string")

        if not sym_tdust:
            raise ValueError("Symbol of dust temperature does not exist")

        if len(spec) != 1:
            raise ValueError("Number of species in thermal desoprtion should be 1.")

        spec = spec[0]
        rate = " * ".join(
            [
                f"opt_thd",
                f"sqrt(2.0*sites*kerg*eb_{spec.alias}/(pi*pi*amu*{spec.massnumber}))",
                f"2.0 * densites",
                f"exp(-eb_{spec.alias}/{sym_tdust})",
            ],
        )
        rate = f"mantabund > 1e-30 ? ({rate}) : 0.0"
        return rate

    def _rate_photon_desorption(
        self,
        spec: list[Species],
        a: float,
        b: float,
        c: float,
        sym_phot: str,
    ) -> str:

        if not isinstance(sym_phot, str):
            raise TypeError("sym_phot should be a string")

        if not sym_phot:
            raise ValueError("Symbol of photon intensity does not exist")

        if len(spec) != 1:
            raise ValueError("Number of species in photon desoprtion should be 1.")

        spec = spec[0]
        rate = " * ".join(
            [
                f"opt_uvd * 4.875e3 * garea",
                f"({sym_phot}) * {spec.photon_yield(default=0.1)} / mant",
            ]
        )

        rate = f"eb_uvd >= {spec.binding_energy} ? ({rate}) : 0.0"
        rate = f"mantabund > 1e-30 ? ({rate}) : 0.0"
        return rate

    def _rate_cosmicray_desorption(
        self,
        spec: list[Species],
        a: float,
        b: float,
        c: float,
        sym_cr: str,
    ) -> str:

        if not isinstance(sym_cr, str):
            raise TypeError("sym_cr should be a string")

        if not sym_cr:
            raise ValueError("Symbol of cosmic ray ionization rate does not exist")

        if len(spec) != 1:
            raise ValueError("Number of species in cosmic-ray desoprtion should be 1.")

        spec = spec[0]
        rate = " * ".join(
            [
                f"opt_crd * 4.0 * pi * crdeseff",
                f"({sym_cr})",
                f"1.64e-4 * garea / mant",
            ]
        )

        rate = f"eb_crd >= {spec.binding_energy} ? ({rate}) : 0.0"
        rate = f"mantabund > 1e-30 ? ({rate}) : 0.0"
        return rate

    def _rate_h2_desorption(
        self,
        spec: list[Species],
        a: float,
        b: float,
        c: float,
        sym_h2form: str,
    ) -> str:

        if not isinstance(sym_h2form, str):
            raise TypeError("sym_h2form should be a string")

        if not sym_h2form:
            raise ValueError("Symbol of H2 formation rate does not exist")

        if len(spec) != 1:
            raise ValueError("Number of species in H2 desoprtion should be 1.")

        spec = spec[0]
        rate = f"opt_h2d * h2deseff * {sym_h2form} / mant"
        rate = f"eb_h2d >= {spec.binding_energy} ? ({rate}) : 0.0"
        rate = f"mantabund > 1e-30 ? ({rate}) : 0.0"

        return rate

    def rateexpr(
        self,
        rtype: ReactionType,
        reactants: list[Species],
        alpha: float = 0.0,
        beta: float = 0.0,
        gamma: float = 0.0,
        sym_tgas: str = "",
        sym_tdust: str = "",
        sym_phot: str = "",
        sym_cr: str = "",
        sym_h2form: str = "",
    ) -> str:

        if rtype == ReactionType.GRAIN_FREEZE:
            return self._rate_depletion(reactants, alpha, beta, gamma, sym_tgas)

        elif rtype == ReactionType.GRAIN_DESORB_THERMAL:
            return self._rate_thermal_desorption(
                reactants, alpha, beta, gamma, sym_tdust
            )

        elif rtype == ReactionType.GRAIN_DESORB_PHOTON:
            return self._rate_photon_desorption(reactants, alpha, beta, gamma, sym_phot)

        elif rtype == ReactionType.GRAIN_DESORB_COSMICRAY:
            return self._rate_cosmicray_desorption(
                reactants, alpha, beta, gamma, sym_cr
            )

        elif rtype == ReactionType.GRAIN_DESORB_H2:
            return self._rate_h2_desorption(reactants, alpha, beta, gamma, sym_h2form)

        else:
            raise ValueError(f"Unknown reaction type in RR07 dust model: {rtype}")
