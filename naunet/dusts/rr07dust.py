from __future__ import annotations
from .dust import Dust
from ..species import Species


class RR07Dust(Dust):

    varis = {
        "Radius": "rG",
        "GrainDensity": "gdens",  # grain density
        "SurfaceSites": "sites",
        "FreezeRatio": "fr",
        "ThermDesorptionOption": "opt_thd",
        "CRDesorptionOption": "opt_crd",
        "H2DesorptionOption": "opt_h2d",
        "UVDesorptionOption": "opt_uvd",
        "CRDesorptionEfficiency": "crdeseff",
        "H2DesorptionEfficiency": "h2deseff",
    }
    user_var = [
        # "double mant = GetMantleDens(y) > 0.0 ? GetMantleDens(y) : 1e-40",
        "double mant = GetMantleDens(y)",
        "double mantabund = mant / nH",
        "double garea = (pi*rG*rG) * gdens",  # total grain cross-section
        "double garea_per_H = garea / nH",
        "double densites = 4.0 * garea * sites",
    ]

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)

    def rate_depletion(
        self, spec: Species, a: float, b: float, c: float, tgas: str
    ) -> str:
        rg = self.varis.get("Radius")
        fr = self.varis.get("FreezeRatio")

        if spec.iselectron:
            rate = f"4.57e4 * {a} * garea * {fr} * ( 1.0 + 16.71e-4/({rg} * {tgas}) )"
        elif spec.charge == 0:
            rate = f"4.57e4 * {a} * sqrt({tgas} / {spec.massnumber}) * garea * {fr}"
        else:
            rate = f"4.57e4 * {a} * sqrt({tgas} / {spec.massnumber}) * garea * {fr} * ( 1.0 + 16.71e-4/({rg} * {tgas}) )"

        return rate

    def rate_desorption(
        self,
        spec: Species,
        a: float,
        b: float,
        c: float,
        tdust: str = "",
        zeta: str = "",
        uvphot: str = "",
        h2form: str = "",
        destype: str = "",
    ) -> str:

        crdeseff = self.varis.get("CRDesorptionEfficiency")
        h2deseff = self.varis.get("H2DesorptionEfficiency")

        if destype == "thermal":

            sites = self.varis.get("SurfaceSites")

            if not tdust:
                raise ValueError("Symbol of dust temperature was not provided.")
            rate = " * ".join(
                [
                    f"opt_thd",
                    f"sqrt(2.0*{sites}*kerg*eb_{spec.alias}/(pi*pi*amu*{spec.massnumber}))",
                    f"2.0 * densites",
                    f"exp(-eb_{spec.alias}/{tdust})",
                ],
            )

        elif destype == "cosmicray":
            if not zeta:
                raise ValueError(
                    "Symbol of cosmic ray ionization rate (in Draine unit) was not provided."
                )
            rate = " * ".join(
                [
                    f"opt_crd * 4.0 * pi * {crdeseff}",
                    f"({zeta})",
                    f"1.64e-4 * garea / mant",
                ]
            )

        elif destype == "photon":
            if not uvphot:
                raise ValueError("Symbol of UV field strength was not provided.")

            rate = " * ".join(
                [
                    f"opt_uvd * 4.875e3 * garea",
                    f"({uvphot}) * {spec.photon_yield(default=0.1)} / mant",
                ]
            )

        elif destype == "h2":
            if not h2form:
                raise ValueError("Symbol of H2 formation rate was not provided.")

            rate = f"opt_h2d * {h2deseff} * {h2form} / mant"

        else:
            raise ValueError(f"Not support desorption type {destype}")

        rate = f"mantabund > 1e-30 ? ({rate}) : 0.0"

        return rate

    def rate_electroncapture(self, tgas) -> str:
        return NotImplemented

    def rate_recombination(self, a: float, b: float, c: float, tgas: str) -> str:
        return NotImplemented

    def rate_surface1(self, *args, **kwargs) -> str:
        return NotImplemented

    def rate_surface2(
        self,
        re1: Species,
        re2: Species,
        a: float,
        b: float,
        c: float,
        tdust: str = "",
        reacdes: bool = False,
    ) -> str:
        return NotImplemented
