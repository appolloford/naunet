from .dust import Dust
from ..species import Species


class UniDust(Dust):

    varis = {
        "Radius": "rG",
        "Barrier": "barr",
        "SurfaceSites": "sites",
        "HOPRatio": "hop",
        "MonoLayers": "nMono",
        "DutyCycle": "duty",
        "CRDesorptionTemperature": "Tcr",
        "BranchRatio": "branch",
    }
    user_var = [
        "double mant = GetMantleDens(y)",
        "double gdens = (y[IDX_GRAIN0I]+y[IDX_GRAINM])",
        "double garea = (4*pi*rG*rG) * gdens",
        "double unisites = sites * (4*pi*rG*rG)",
        "double densites = sites * (4*pi*rG*rG) * gdens",
        "double freq = sqrt((2.0*sites*kerg)/((pi*pi)*amu))",
        "double quan = -2.0*(barr/hbar) * sqrt(2.0*amu*kerg)",
        "double layers = mant/(nMono*densites)",
        "double cov = (mant == 0.0) ? 0.0 : fmin(layers/mant, 1.0/mant)",
    ]

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)

    def rate_depletion(self, a: float, b: float, c: float, tgas: str) -> str:
        rg = self.varis.get("Radius")
        rate = f"{a} * pi * pow({rg}, 2.0) * gdens * sqrt(8.0 * kerg * {tgas}/ (pi*amu*{c}))"
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
        destype: str = "",
    ) -> str:

        sites = self.varis.get("SurfaceSites")
        nmono = self.varis.get("MonoLayers")
        duty = self.varis.get("DutyCycle")
        Tcr = self.varis.get("CRDesorptionTemperature")

        if destype == "thermal":
            if not tdust:
                raise ValueError("Symbol of dust temperature was not provided.")
            rate = " * ".join(
                [
                    f"sqrt(2.0*{sites}*kerg*eb_{spec.alias}/(pi*pi*amu*{c}))",
                    f"{nmono}",
                    "densites",
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
                    f"({zeta})",
                    f"{duty}",
                    f"sqrt(2.0*{sites}*kerg*eb_{spec.alias}/(pi*pi*amu*{c}))",
                    f"{nmono}",
                    "densites",
                    f"exp(-eb_{spec.alias}/{Tcr})",
                ]
            )
        elif destype == "photon":
            if not uvphot:
                raise ValueError("Symbol of UV field strength (G0) was not provided.")
            rate = f"({uvphot}) * {spec.photon_yield} * {nmono} * garea"
        else:
            raise ValueError(f"Not support desorption type {destype}")

        rate = " * ".join([rate, "cov"])

        return rate

    def rate_electroncapture(self, tgas) -> str:
        rg = self.varis.get("Radius")
        return f"pi * {rg} * {rg} * sqrt(8.0*kerg*{tgas}/pi/amu/meu)"

    def rate_recombination(self, a: float, b: float, c: float, tgas: str) -> str:
        rg = self.varis.get("Radius")
        rate = " * ".join(
            [
                f"{a} * pi * pow({rg}, 2.0) * gdens",
                f"sqrt(8.0*kerg*{tgas}/(pi*amu*{c}))",
                f"(1.0 + pow(echarge, 2.0)/{rg}/kerg/{tgas})",
                f"(1.0 + sqrt(2.0*pow(echarge, 2.0)/({rg}*kerg*{tgas}+2.0*pow(echarge, 2.0))))",
            ]
        )
        return rate

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
        hop = self.varis.get("HOPRatio")
        nmono = self.varis.get("MonoLayers")
        eb1 = re1.binding_energy
        eb2 = re2.binding_energy
        nmass1 = re1.massnumber
        nmass2 = re2.massnumber

        afreq = f"freq * sqrt({eb1}/{nmass1})"
        adiff = f"{afreq} * exp(-{eb1}*{hop}/{tdust})/unisites"
        aquan = f"{afreq} * exp(quan * sqrt({hop}*{nmass1}*{eb1})) / unisites"

        bfreq = f"freq * sqrt({eb2}/{nmass2})"
        bdiff = f"{bfreq} * exp(-{eb2}*{hop}/{tdust})/unisites"
        bquan = f"{bfreq} * exp(quan * sqrt({hop}*{nmass2}*{eb2})) / unisites"

        kappa = f"exp(-{a}/{tdust})"
        kquan = f"exp(quan * sqrt((({nmass1}*{nmass2})/({nmass1}+{nmass2}))*{a}))"

        if re1.name in ["GH", "GH2"] and re2.name in ["GH", "GH2"]:
            rate = " * ".join(
                [
                    f"fmax({kappa}, {kquan})",
                    f"(fmax({adiff}, {aquan})+fmax({bdiff}, {bquan}))",
                    f"pow(({nmono}*densites), 2.0) / gdens",
                ]
            )
        elif re1.name in ["GH", "GH2"]:
            rate = " * ".join(
                [
                    f"fmax({kappa}, {kquan})",
                    f"(fmax({adiff}, {aquan})+{bdiff})",
                    f"pow(({nmono}*densites), 2.0) / gdens",
                ]
            )
        elif re2.name in ["GH", "GH2"]:
            rate = " * ".join(
                [
                    f"fmax({kappa}, {kquan})",
                    f"({adiff}+fmax({bdiff}, {bquan}))",
                    f"pow(({nmono}*densites), 2.0) / gdens",
                ]
            )
        else:
            rate = " * ".join(
                [
                    f"{kappa} * ({adiff}+{bdiff})",
                    f"pow(({nmono}*densites), 2.0) / gdens",
                ]
            )

        rate = " * ".join([rate, "cov", "cov"])

        if reacdes:
            branch = self.varis.get("BranchRatio")
            rate = f"{branch} * {rate}"

        return rate
