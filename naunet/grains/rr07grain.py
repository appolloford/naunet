from __future__ import annotations
from typing import TYPE_CHECKING
from .grain import Grain
from ..species import Species
from ..reactions.reactiontype import ReactionType

if TYPE_CHECKING:
    from ..reactions.reaction import Reaction


class RR07Grain(Grain):
    """
    The dust grain model from Roberts et al. (2007). Follow the implementation in
    UCLCHEMv1.3.

    The class inherit class `Grain`. For the general description, please refer to the
    page of `Grain`.

    Attributes:
        model (str): The name of the model = "rr07"
        varis (dict): The parameters exist in the model
        locvars (list): The local variables in the model
    """

    model = "rr07"

    varis = {
        "rG": 1e-5,  # grain radius
        "gdens": 7.6394373e-13,  # grain density
        "sites": 1e15,  # surface sites
        "fr": 1.0,  # freeze ratio
        "opt_crd": 1.0,  # cosmic ray induced desorption option
        "opt_h2d": 1.0,  # H2 formation induced desorption option
        "opt_uvd": 1.0,  # UV desorption option
        "eb_h2d": 1.21e3,  # maxmium binding energy which H2 desorption can desorb
        "eb_crd": 1.21e3,  # maxmium binding energy which CR desorption can desorb
        "eb_uvd": 1.0e4,  # maxmium binding energy which UV desorption can desorb
        "crdeseff": 1e5,  # cosmic ray desorption efficiency
        "h2deseff": 1.0e-2,  # H2 desorption efficiency
        "uvcreff": 1.0e-3,  # UVCREFF is ratio of CR induced UV to ISRF UV
    }
    locvars = [
        # "double mant = GetMantleDens(y) > 0.0 ? GetMantleDens(y) : 1e-40",
        "double mant = GetMantleDens(y)",
        "double mantabund = mant / nH",
        "double gxsec = (pi*rG*rG) * gdens",  # total grain cross-section
        "double garea = 4.0 * (pi*rG*rG) * gdens",
        "double gxsec_per_H = gxsec / nH",
        "double unisites = sites * (4*pi*rG*rG)",
        "double densites = garea * sites",
    ]

    def rate_depletion(self, reac: Reaction) -> str:

        super().rate_depletion(reac)

        spec = reac.reactants[0]
        a = reac.alpha
        tgas = reac.symbols.temperature.symbol

        if spec.is_electron:
            rate = " * ".join(
                [
                    f"4.57e4 * {a} * gxsec * fr",
                    f"( 1.0 + 16.71e-4/(rG * {tgas}) )",
                ]
            )
        elif spec.charge == 0:
            rate = " * ".join(
                [
                    f"4.57e4 * {a} * gxsec * fr",
                    f"sqrt({tgas} / {spec.massnumber})",
                ]
            )
        else:
            rate = " * ".join(
                [
                    f"4.57e4 * {a} * gxsec * fr",
                    f"sqrt({tgas} / {spec.massnumber})",
                    f"( 1.0 + 16.71e-4/(rG * {tgas}) )",
                ]
            )

        return rate

    def rate_photon_desorption(self, reac: Reaction) -> str:

        super().rate_photon_desorption(reac)

        crrate = reac.symbols.cosmic_ray_ionization_rate.symbol
        radfield = reac.symbols.radiation_field
        av = reac.symbols.visual_extinction

        if not self.sym_crrate or not self.sym_av:
            raise ValueError(
                "Cosmic-ray ionization rate symbol or visual"
                "extinction symbol is not set properly"
            )

        sym_phot = f"({crrate} + ({radfield} / uvcreff) * exp(-1.8*{av}) )"

        spec = reac.reactants[0]
        rate = " * ".join(
            [
                f"opt_uvd * 4.875e3 * gxsec",
                f"({sym_phot}) * {spec.photon_yield(default=0.1)} / mant",
            ]
        )

        rate = f"eb_uvd >= {spec.binding_energy} ? ({rate}) : 0.0"
        rate = f"mantabund > 1e-30 ? ({rate}) : 0.0"
        return rate

    def rate_cosmicray_desorption(self, reac: Reaction) -> str:

        super().rate_cosmicray_desorption(reac)

        crrate = reac.symbols.cosmic_ray_ionization_rate.symbol

        spec = reac.reactants[0]
        rate = " * ".join(
            [
                f"opt_crd * 4.0 * pi * crdeseff",
                f"({crrate})",
                f"1.64e-4 * gxsec / mant",
            ]
        )

        rate = f"eb_crd >= {spec.binding_energy} ? ({rate}) : 0.0"
        rate = f"mantabund > 1e-30 ? ({rate}) : 0.0"
        return rate

    def rate_h2_desorption(self, reac: Reaction) -> str:

        super().rate_h2_desorption(reac)

        h2form = reac.symbols.H2_formation_rate.symbol

        spec = reac.reactants[0]
        rate = f"opt_h2d * h2deseff * {h2form} / mant"
        rate = f"eb_h2d >= {spec.binding_energy} ? ({rate}) : 0.0"
        rate = f"mantabund > 1e-30 ? ({rate}) : 0.0"

        return rate


class RR07XGrain(RR07Grain):

    model = "rr07x"

    varis = {
        **RR07Grain.varis,
        "opt_thd": 1.0,  # thermal desorption option
    }

    def rate_thermal_desorption(self, reac: Reaction) -> str:

        super().rate_thermal_desorption(reac)

        tdust = reac.symbols.dust_temperature.symbol

        spec = reac.reactants[0]
        rate = " * ".join(
            [
                f"opt_thd",
                f"sqrt(2.0*sites*kerg*eb_{spec.alias}/(pi*pi*amu*{spec.massnumber}))",
                f"2.0 * densites",
                f"exp(-eb_{spec.alias}/{tdust})",
            ],
        )
        rate = f"mantabund > 1e-30 ? ({rate}) : 0.0"
        return rate
