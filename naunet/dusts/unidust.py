from .dust import Dust
from ..species import Species


class UniDust(Dust):

    varis = {
        "rG": 1.0e-5,  # grain radius
        "barr": 1.5e-8,  # barrier
        "sites": 1e15,  # surface sites
        "hop": 0.3,  # hop ration
        "nMono": 2.0,  # number of monolayers
        "duty": 3.16e-19,  # duty cycle
        "Tcr": 70.0,  # cosmic ray induced desorption temperature
        "branch": 1e-2,  # branch ration
    }
    locvars = [
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
        rate = f"{a} * pi * rG * rG * gdens * sqrt(8.0 * kerg * {tgas}/ (pi*amu*{c}))"
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

        if destype == "thermal":
            if not tdust:
                raise ValueError("Symbol of dust temperature was not provided.")
            rate = " * ".join(
                [
                    f"sqrt(2.0*sites*kerg*eb_{spec.alias}/(pi*pi*amu*{c}))",
                    f"nMono",
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
                    f"duty",
                    f"sqrt(2.0*sites*kerg*eb_{spec.alias}/(pi*pi*amu*{c}))",
                    f"nMono",
                    "densites",
                    f"exp(-eb_{spec.alias}/Tcr)",
                ]
            )
        elif destype == "photon":
            if not uvphot:
                raise ValueError("Symbol of UV field strength (G0) was not provided.")
            rate = f"({uvphot}) * {spec.photon_yield()} * nMono * garea"
        else:
            raise ValueError(f"Not support desorption type {destype}")

        rate = " * ".join([rate, "cov"])

        return rate

    def rate_electroncapture(self, tgas) -> str:
        return f"pi * rG * rG * sqrt(8.0*kerg*{tgas}/pi/amu/meu)"

    def rate_recombination(self, a: float, b: float, c: float, tgas: str) -> str:
        rate = " * ".join(
            [
                f"{a} * pi * rG * rG * gdens",
                f"sqrt(8.0*kerg*{tgas}/(pi*amu*{c}))",
                f"(1.0 + pow(echarge, 2.0)/rG/kerg/{tgas})",
                f"(1.0 + sqrt(2.0*pow(echarge, 2.0)/(rG*kerg*{tgas}+2.0*pow(echarge, 2.0))))",
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

        eb1 = re1.binding_energy
        eb2 = re2.binding_energy
        nmass1 = re1.massnumber
        nmass2 = re2.massnumber

        afreq = f"freq * sqrt({eb1}/{nmass1})"
        adiff = f"{afreq} * exp(-{eb1}*hop/{tdust})/unisites"
        aquan = f"{afreq} * exp(quan * sqrt(hop*{nmass1}*{eb1})) / unisites"

        bfreq = f"freq * sqrt({eb2}/{nmass2})"
        bdiff = f"{bfreq} * exp(-{eb2}*hop/{tdust})/unisites"
        bquan = f"{bfreq} * exp(quan * sqrt(hop*{nmass2}*{eb2})) / unisites"

        kappa = f"exp(-{a}/{tdust})"
        kquan = f"exp(quan * sqrt((({nmass1}*{nmass2})/({nmass1}+{nmass2}))*{a}))"

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

        if reacdes:

            rate = f"branch * {rate}"

        return rate
