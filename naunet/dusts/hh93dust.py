from __future__ import annotations
from .dust import Dust
from ..species import Species
from ..reactions.reactiontype import ReactionType


class HH93Dust(Dust):

    varis = {
        "rG": 1.0e-5,  # grain radius
        "barr": 1.5e-8,  # barrier
        "sites": 1e15,  # surface sites
        "hop": 0.3,  # hop ratio
        "nMono": 2.0,  # number of monolayers
        "duty": 3.16e-19,  # duty cycle
        "Tcr": 70.0,  # cosmic ray induced desorption temperature
        "branch": 1e-2,  # branch ratio
        "opt_frz": 1.0,
        "opt_thd": 1.0,
        "opt_uvd": 1.0,
        "opt_crd": 1.0,
        "opt_rcd": 1.0,  # reactive desorption option
    }
    locvars = [
        "double mant = GetMantleDens(y)",
        "double garea = (4*pi*rG*rG) * gdens",
        "double unisites = sites * (4*pi*rG*rG)",
        "double densites = sites * (4*pi*rG*rG) * gdens",
        "double freq = sqrt((2.0*sites*kerg)/((pi*pi)*amu))",
        "double quan = -2.0*(barr/hbar) * sqrt(2.0*amu*kerg)",
        "double layers = mant/(nMono*densites)",
        "double cov = (mant == 0.0) ? 0.0 : fmin(layers/mant, 1.0/mant)",
    ]

    def __init__(self, species: list[str] | list[Species] = None) -> None:
        super().__init__(model="hh93", species=species)

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

        rate = " * ".join(
            [
                f"opt_frz * {a} * pi * rG * rG * gdens",
                f"sqrt(8.0 * kerg * {sym_tgas}/ (pi*amu*{c}))",
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
                f"cov",
                f"opt_thd * nMono * densites",
                f"sqrt(2.0*sites*kerg*eb_{spec.alias}/(pi*pi*amu*{c}))",
                f"exp(-eb_{spec.alias}/({sym_tdust}))",
            ],
        )
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
        rate = f"cov * ({sym_phot}) * {spec.photon_yield()} * nMono * garea"
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
                f"cov",
                f"duty * nMono * densites",
                f"({sym_cr})",
                f"sqrt(2.0*sites*kerg*eb_{spec.alias}/(pi*pi*amu*{c}))",
                f"exp(-eb_{spec.alias}/Tcr)",
            ]
        )
        return rate

    def _rate_electron_capture(self, sym_tgas: str) -> str:

        if not isinstance(sym_tgas, str):
            raise TypeError("sym_tgas should be a string")

        if not sym_tgas:
            raise ValueError("Symbol of gas temperature does not exist")

        rate = f"pi * rG * rG * sqrt(8.0*kerg*({sym_tgas})/pi/amu/meu)"
        return rate

    def _rate_recombination(
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

        rate = " * ".join(
            [
                f"{a} * pi * rG * rG * gdens",
                f"sqrt(8.0*kerg*{sym_tgas}/(pi*amu*{c}))",
                f"(1.0 + pow(echarge, 2.0)/rG/kerg/{sym_tgas})",
                f"(1.0 + sqrt(2.0*pow(echarge, 2.0)/(rG*kerg*{sym_tgas}+2.0*pow(echarge, 2.0))))",
            ]
        )
        return rate

    def _rate_surface_twobody(
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

        if len(spec) != 2:
            raise ValueError(
                "Number of species in two-body surface reaction should be 2."
            )

        re1, re2 = spec
        eb1, nmass1, eb2, nmass2 = re1.eb, re1.A, re2.eb, re2.A

        afreq = f"freq * sqrt({eb1}/{nmass1})"
        adiff = f"{afreq} * exp(-{eb1}*hop/{sym_tdust})/unisites"
        aquan = f"{afreq} * exp(quan * sqrt(hop*{nmass1}*{eb1})) / unisites"

        bfreq = f"freq * sqrt({eb2}/{nmass2})"
        bdiff = f"{bfreq} * exp(-{eb2}*hop/{sym_tdust})/unisites"
        bquan = f"{bfreq} * exp(quan * sqrt(hop*{nmass2}*{eb2})) / unisites"

        kappa = f"exp(-{a}/{sym_tdust})"
        kquan = f"exp(quan * sqrt((({nmass1}*{nmass2})/({nmass1}+{nmass2}))*{a}))"

        rate = ""
        if re1.name in ["GH", "GH2"] and re2.name in ["GH", "GH2"]:
            rate = " * ".join(
                [
                    f"fmax({kappa}, {kquan})",
                    f"(fmax({adiff}, {aquan})+fmax({bdiff}, {bquan}))",
                    f"pow((nMono*densites), 2.0) / gdens",
                ]
            )
        elif re1.name in ["GH", "GH2"]:
            rate = " * ".join(
                [
                    f"fmax({kappa}, {kquan})",
                    f"(fmax({adiff}, {aquan})+{bdiff})",
                    f"pow((nMono*densites), 2.0) / gdens",
                ]
            )
        elif re2.name in ["GH", "GH2"]:
            rate = " * ".join(
                [
                    f"fmax({kappa}, {kquan})",
                    f"({adiff}+fmax({bdiff}, {bquan}))",
                    f"pow((nMono*densites), 2.0) / gdens",
                ]
            )
        else:
            rate = " * ".join(
                [
                    f"{kappa} * ({adiff}+{bdiff})",
                    f"pow((nMono*densites), 2.0) / gdens",
                ]
            )

        rate = " * ".join([rate, "cov", "cov"])

        return rate

    def _rate_reactive_desorption(
        self,
        spec: list[Species],
        a: float,
        b: float,
        c: float,
        sym_tdust: str = "",
    ) -> str:

        rate = self._rate_surface_twobody(spec, a, b, c, sym_tdust)
        rate = f"branch * {rate}"

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
    ) -> str:

        if rtype == ReactionType.GAS_LEEDS_RECOM:
            return self._rate_recombination(reactants, alpha, beta, gamma, sym_tgas)

        elif rtype == ReactionType.GRAIN_FREEZE:
            return self._rate_depletion(reactants, alpha, beta, gamma, sym_tgas)

        elif rtype == ReactionType.GRAIN_DESORPT_THERMAL:
            return self._rate_thermal_desorption(
                reactants, alpha, beta, gamma, sym_tdust
            )

        elif rtype == ReactionType.GRAIN_DESORPT_PHOTON:
            return self._rate_photon_desorption(reactants, alpha, beta, gamma, sym_phot)

        elif rtype == ReactionType.GRAIN_DESORPT_COSMICRAY:
            return self._rate_cosmicray_desorption(
                reactants, alpha, beta, gamma, sym_cr
            )

        elif rtype == ReactionType.SURFACE_TWOBODY:
            return self._rate_surface_twobody(reactants, alpha, beta, gamma, sym_tdust)

        elif rtype == ReactionType.GRAIN_DESORPT_REACTIVE:
            return self._rate_reactive_desorption(
                reactants, alpha, beta, gamma, sym_tdust
            )

        elif rtype == ReactionType.GAS_LEEDS_ECAPTURE:
            return self._rate_electron_capture(sym_tgas)

        else:
            raise ValueError(f"Unknown reaction type in HH93 dust model: {rtype}")
