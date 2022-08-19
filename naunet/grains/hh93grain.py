from __future__ import annotations
from typing import TYPE_CHECKING
from .grain import Grain
from ..species import Species
from ..reactions.reactiontype import ReactionType


if TYPE_CHECKING:
    from ..reactions.reaction import Reaction


class HH93Grain(Grain):

    model = "hh93"

    consts = {
        "habing": 1e8,
        "crphot": 1e4,
        "hbar": 1.054571726e-27,
    }

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
        super().__init__(species=species)

    def rate_depletion(self, reac: Reaction) -> str:

        super().rate_depletion(reac)

        spec = reac.reactants[0]
        a = reac.alpha

        rate = " * ".join(
            [
                f"opt_frz * {a} * pi * rG * rG * gdens",
                f"sqrt(8.0 * kerg * {self.sym_tgas}/ (pi*amu*{spec.A}))",
            ]
        )
        return rate

    def rate_thermal_desorption(self, reac: Reaction) -> str:

        super().rate_thermal_desorption(reac)

        spec = reac.reactants[0]
        rate = " * ".join(
            [
                f"opt_thd * cov",
                f"nMono * densites",
                f"sqrt(2.0*sites*kerg*eb_{spec.alias}/(pi*pi*amu*{spec.A}))",
                f"exp(-eb_{spec.alias}/({self.sym_tdust}))",
            ],
        )
        return rate

    def rate_photon_desorption(self, reac: Reaction) -> str:

        super().rate_photon_desorption(reac)

        if not self.sym_crrate or not self.sym_av:
            raise ValueError(
                "Cosmic-ray ionization rate symbol or visual"
                "extinction symbol is not set properly"
            )

        sym_phot = f"{self.sym_radfield}*habing*exp(-{self.sym_av}*3.02) + crphot * ({self.sym_crrate})"

        spec = reac.reactants[0]
        rate = f"opt_uvd * cov * ({sym_phot}) * {spec.photon_yield()} * nMono * garea"
        return rate

    def rate_cosmicray_desorption(self, reac: Reaction) -> str:

        super().rate_cosmicray_desorption(reac)

        spec = reac.reactants[0]
        rate = " * ".join(
            [
                f"opt_crd * cov",
                f"duty * nMono * densites",
                f"({self.sym_crrate})",
                f"sqrt(2.0*sites*kerg*eb_{spec.alias}/(pi*pi*amu*{spec.A}))",
                f"exp(-eb_{spec.alias}/Tcr)",
            ]
        )
        return rate

    def rate_electron_capture(self, reac: Reaction) -> str:

        super().rate_electron_capture(reac)

        rate = f"pi * rG * rG * sqrt(8.0*kerg*({self.sym_tgas})/pi/amu/meu)"
        return rate

    def rate_recombination(self, reac: Reaction) -> str:

        super().rate_recombination(reac)

        [spec] = [s for s in reac.reactants if not s.is_grain]
        a = reac.alpha

        rate = " * ".join(
            [
                f"{a} * pi * rG * rG * gdens",
                f"sqrt(8.0*kerg*{self.sym_tgas}/(pi*amu*{spec.A}))",
                f"(1.0 + pow(echarge, 2.0)/rG/kerg/{self.sym_tgas})",
                f"(1.0 + sqrt(2.0*pow(echarge, 2.0)/(rG*kerg*{self.sym_tgas}+2.0*pow(echarge, 2.0))))",
            ]
        )
        return rate

    def _rate_surface(self, reac: Reaction) -> str:

        re1, re2 = reac.reactants
        eb1, nmass1, eb2, nmass2 = re1.eb, re1.A, re2.eb, re2.A
        a = reac.alpha

        afreq = f"freq * sqrt({eb1}/{nmass1})"
        adiff = f"{afreq} * exp(-{eb1}*hop/{self.sym_tdust})/unisites"
        aquan = f"{afreq} * exp(quan * sqrt(hop*{nmass1}*{eb1})) / unisites"

        bfreq = f"freq * sqrt({eb2}/{nmass2})"
        bdiff = f"{bfreq} * exp(-{eb2}*hop/{self.sym_tdust})/unisites"
        bquan = f"{bfreq} * exp(quan * sqrt(hop*{nmass2}*{eb2})) / unisites"

        kappa = f"exp(-{a}/{self.sym_tdust})"
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

    def rate_surface_twobody(self, reac: Reaction) -> str:

        super().rate_surface_twobody(reac)

        rate = self._rate_surface(reac)
        return rate

    def rate_reactive_desorption(self, reac: Reaction) -> str:

        super().rate_reactive_desorption(reac)

        rate = self._rate_surface(reac)
        rate = f"opt_rcd * branch * {rate}"

        return rate
