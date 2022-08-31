from __future__ import annotations
from typing import TYPE_CHECKING
from .grain import Grain
from ..component import VariableType as vt
from ..species import Species
from ..reactiontype import ReactionType

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

    def __init__(self, species: list[Species] = None, group: int = 0) -> None:
        super().__init__(species, group)

        group = group or ""
        self.register(
            "grain_density",
            (f"gdens{group}", 7.6394373e-13, vt.param),
            force_overwrite=True,
        )
        self.register(
            "grain_surface_sites_density",
            (f"sites{group}", 1.5e15, vt.param),
        )
        self.register("freeze_option", (f"fr{group}", 1.0, vt.param))
        self.register(
            "cosmic_ray_desorption_option",
            (f"opt_crd{group}", 1.0, vt.param),
        )
        self.register("photon_desorption_option", (f"opt_uvd{group}", 1.0, vt.param))
        self.register("H2_desorption_option", (f"opt_h2d{group}", 1.0, vt.param))
        self.register(
            "max_cosmic_ray_desorption_binding_energy",
            (f"eb_crd{group}", 1.21e3, vt.param),
        )
        self.register(
            "max_photon_desorption_binding_energy",
            (f"eb_uvd{group}", 1.0e4, vt.param),
        )
        self.register(
            "max_H2_desorption_binding_energy",
            (f"eb_h2d{group}", 1.21e3, vt.param),
        )
        self.register(
            "cosmic_ray_desorption_efficiency",
            (f"crdeseff{group}", 1.0e5, vt.param),
        )
        self.register(
            "cosmic_ray_induce_photon_efficiency",
            (f"uvcreff{group}", 1.0e-3, vt.param),
        )
        self.register(
            "H2_desorption_efficiency",
            (f"h2deseff{group}", 1.0e-2, vt.param),
        )

        # TODO: get mantle density by group
        self.register(
            "mantle_number_density",
            (f"mant{group}", "GetMantleDens(y)", vt.derived),
        )
        self.register(
            "mantle_number_density_per_H",
            (f"mantabund{group}", f"mant{group} / nH", vt.derived),
        )
        self.register(
            "grain_cross_section",
            (f"gxsec{group}", f"(pi*rG{group}*rG{group}) * gdens", vt.derived),
        )
        self.register(
            "grain_surface_area",
            (f"garea{group}", f"4.0 * gxsec{group}", vt.derived),
        )
        self.register(
            "unit_surface_sites",
            (
                f"unisites{group}",
                f"sites{group} * (4*pi*rG{group}*rG{group})",
                vt.derived,
            ),
        )
        self.register(
            "surface_sites_density",
            (f"densites{group}", f"garea{group} * sites{group}", vt.derived),
        )

    def rate_depletion(self, reac: Reaction) -> str:

        super().rate_depletion(reac)

        spec = reac.reactants[0]
        a = reac.alpha
        tgas = reac.symbols.temperature.symbol

        r = self.symbols.grain_radius.symbol
        gxsec = self.symbols.grain_cross_section.symbol
        fr = self.symbols.freeze_option.symbol

        if spec.is_electron:
            rate = " * ".join(
                [
                    f"4.57e4 * {a} * {gxsec} * {fr}",
                    f"( 1.0 + 16.71e-4/({r} * {tgas}) )",
                ]
            )
        elif spec.charge == 0:
            rate = " * ".join(
                [
                    f"4.57e4 * {a} * {gxsec} * {fr}",
                    f"sqrt({tgas} / {spec.massnumber})",
                ]
            )
        else:
            rate = " * ".join(
                [
                    f"4.57e4 * {a} * {gxsec} * {fr}",
                    f"sqrt({tgas} / {spec.massnumber})",
                    f"( 1.0 + 16.71e-4/({r} * {tgas}) )",
                ]
            )

        return rate

    def rate_photon_desorption(self, reac: Reaction) -> str:

        super().rate_photon_desorption(reac)

        crrate = reac.symbols.cosmic_ray_ionization_rate.symbol
        zism = reac.symbols.ism_cosmic_ray_ionization_rate.symbol
        radfield = reac.symbols.radiation_field
        av = reac.symbols.visual_extinction

        mant = self.symbols.mantle_number_density.symbol
        mantabund = self.symbols.mantle_number_density_per_H.symbol
        gxsec = self.symbols.grain_cross_section.symbol
        opt_uvd = self.symbols.photon_desorption_option.symbol
        eb_uvd = self.symbols.max_photon_desorption_binding_energy.symbol
        uvcreff = self.symbols.cosmic_ray_induce_photon_efficiency.symbol

        sym_phot = f"(({crrate} / {zism}) + ({radfield} / {uvcreff}) * exp(-1.8*{av}) )"

        spec = reac.reactants[0]
        rate = " * ".join(
            [
                f"{opt_uvd} * 4.875e3 * {gxsec}",
                f"({sym_phot}) * {spec.photon_yield(default=0.1)} / {mant}",
            ]
        )

        rate = f"{eb_uvd} >= {spec.binding_energy} ? ({rate}) : 0.0"
        rate = f"{mantabund} > 1e-30 ? ({rate}) : 0.0"
        return rate

    def rate_cosmicray_desorption(self, reac: Reaction) -> str:

        super().rate_cosmicray_desorption(reac)

        crrate = reac.symbols.cosmic_ray_ionization_rate.symbol
        zism = reac.symbols.ism_cosmic_ray_ionization_rate.symbol

        mant = self.symbols.mantle_number_density.symbol
        mantabund = self.symbols.mantle_number_density_per_H.symbol
        gxsec = self.symbols.grain_cross_section.symbol
        opt_crd = self.symbols.cosmic_ray_desorption_option.symbol
        eb_crd = self.symbols.max_cosmic_ray_desorption_binding_energy.symbol
        crdeseff = self.symbols.cosmic_ray_desorption_efficiency.symbol

        spec = reac.reactants[0]
        rate = " * ".join(
            [
                f"{opt_crd} * 4.0 * pi * {crdeseff}",
                f"({crrate} / {zism})",
                f"1.64e-4 * {gxsec} / {mant}",
            ]
        )

        rate = f"{eb_crd} >= {spec.binding_energy} ? ({rate}) : 0.0"
        rate = f"{mantabund} > 1e-30 ? ({rate}) : 0.0"
        return rate

    def rate_h2_desorption(self, reac: Reaction) -> str:

        super().rate_h2_desorption(reac)

        h2form = reac.symbols.H2_formation_rate.symbol

        mant = self.symbols.mantle_number_density.symbol
        mantabund = self.symbols.mantle_number_density_per_H.symbol
        opt_h2d = self.symbols.H2_desorption_option.symbol
        eb_h2d = self.symbols.max_H2_desorption_binding_energy.symbol
        h2deseff = self.symbols.H2_desorption_efficiency.symbol

        spec = reac.reactants[0]
        rate = f"{opt_h2d} * {h2deseff} * {h2form} * y[IDX_HI] / {mant}"
        rate = f"{eb_h2d} >= {spec.binding_energy} ? ({rate}) : 0.0"
        rate = f"{mantabund} > 1e-30 ? ({rate}) : 0.0"

        return rate


class RR07XGrain(RR07Grain):

    model = "rr07x"

    def __init__(self, species: list[Species] = None, group: int = 0) -> None:
        super().__init__(species, group)

        self.register("thermal_desorption_option", ("opt_thd", 1.0, vt.param))

    def rate_thermal_desorption(self, reac: Reaction) -> str:

        super().rate_thermal_desorption(reac)

        tdust = reac.symbols.dust_temperature.symbol

        mantabund = self.symbols.mantle_number_density_per_H.symbol
        sites = self.symbols.grain_surface_sites_density.symbol
        densites = self.symbols.surface_sites_density.symbol
        opt_thd = self.symbols.thermal_desorption_option.symbol

        spec = reac.reactants[0]
        rate = " * ".join(
            [
                f"{opt_thd}",
                f"sqrt(2.0*{sites}*kerg*eb_{spec.alias}/(pi*pi*amu*{spec.massnumber}))",
                f"2.0 * {densites}",
                f"exp(-eb_{spec.alias}/{tdust})",
            ],
        )
        rate = f"{mantabund} > 1e-30 ? ({rate}) : 0.0"
        return rate
