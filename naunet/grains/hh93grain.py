from __future__ import annotations
from typing import TYPE_CHECKING
from .grain import Grain
from ..component import VariableType as vt
from ..species import Species
from ..reactiontype import ReactionType


if TYPE_CHECKING:
    from ..reactions.reaction import Reaction


class HH93Grain(Grain):

    model = "hh93"

    def __init__(self, species: list[Species] = None, group: int = 0) -> None:
        super().__init__(species, group)

        group = group or ""
        self.register("habing_field_photon_number", ("habing", 1e8, vt.constant))
        self.register("cosmic_ray_induced_photon_number", ("crphot", 1e4, vt.constant))
        self.register(
            "grain_surface_sites_density",
            (f"sites{group}", 1.5e15, vt.param),
        )
        self.register("grain_sites_barrier_height", (f"barr{group}", 1.5e-8, vt.param))
        self.register("surface_hopping_ratio", (f"hop{group}", 0.3, vt.param))
        self.register("active_monolayer_number", (f"nMono{group}", 2.0, vt.param))
        self.register("freeze_option", (f"opt_frz{group}", 1.0, vt.param))
        self.register("thermal_desorption_option", (f"opt_thd{group}", 1.0, vt.param))
        self.register(
            "cosmic_ray_desorption_option",
            (f"opt_crd{group}", 1.0, vt.param),
        )
        self.register(
            "cosmic_ray_desorption_duty_cycle", (f"duty{group}", 3.16e-19, vt.param)
        )
        self.register(
            "max_cosmic_ray_desorption_temperature", (f"Tcr{group}", 70.0, vt.param)
        )
        self.register("photon_desorption_option", (f"opt_uvd{group}", 1.0, vt.param))
        self.register("reactive_desorption_option", (f"opt_rcd{group}", 1.0, vt.param))
        self.register(
            "reactive_desorption_branching_ratio",
            (f"branch{group}", 1e-2, vt.param),
        )

        # TODO: get mantle density by group
        self.register(
            "mantle_number_density",
            (f"mant{group}", "GetMantleDens(y)", vt.derived),
        )
        self.register(
            "grain_surface_area",
            (f"garea{group}", f"(4.0*pi*rG{group}*rG{group}) * gdens", vt.derived),
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
        self.register(
            "classical_diffusion_rate_factor",
            (
                f"freq{group}",
                f"sqrt((2.0*sites{group}*kerg)/((pi*pi)*amu))",
                vt.derived,
            ),
        )
        self.register(
            "quantum_diffusion_rate_factor",
            (
                f"quan{group}",
                f"-2.0*(barr{group}/hbar) * sqrt(2.0*amu*kerg)",
                vt.derived,
            ),
        )
        self.register(
            "monolayer_number",
            (
                f"layers{group}",
                f"mant{group}/(nMono{group}*densites{group})",
                vt.derived,
            ),
        )
        self.register(
            "coverage",
            (
                f"cov{group}",
                f"(mant{group} == 0.0) ? 0.0 : fmin(layers{group}/mant{group}, 1.0/mant{group})",
                vt.derived,
            ),
        )

    def rate_depletion(self, reac: Reaction) -> str:

        super().rate_depletion(reac)

        spec = reac.reactants[0]
        a = reac.alpha
        tgas = reac.symbols.temperature.symbol

        r = self.symbols.grain_radius.symbol
        gdens = self.symbols.grain_density.symbol
        opt_frz = self.symbols.freeze_option.symbol

        rate = " * ".join(
            [
                f"{opt_frz} * {a} * pi * {r} * {r} * {gdens}",
                f"sqrt(8.0 * kerg * {tgas}/ (pi*amu*{spec.A}))",
            ]
        )
        return rate

    def rate_thermal_desorption(self, reac: Reaction) -> str:

        super().rate_thermal_desorption(reac)

        tdust = reac.symbols.dust_temperature.symbol

        opt_thd = self.symbols.thermal_desorption_option.symbol
        cov = self.symbols.coverage.symbol
        nMono = self.symbols.active_monolayer_number.symbol
        sites = self.symbols.grain_surface_sites_density.symbol
        densites = self.symbols.surface_sites_density.symbol

        spec = reac.reactants[0]
        rate = " * ".join(
            [
                f"{opt_thd} * {cov}",
                f"{nMono} * {densites}",
                f"sqrt(2.0*{sites}*kerg*eb_{spec.alias}/(pi*pi*amu*{spec.A}))",
                f"exp(-eb_{spec.alias}/({tdust}))",
            ],
        )
        return rate

    def rate_photon_desorption(self, reac: Reaction) -> str:

        super().rate_photon_desorption(reac)

        crrate = reac.symbols.cosmic_ray_ionization_rate.symbol
        zism = reac.symbols.ism_cosmic_ray_ionization_rate.symbol
        radfield = reac.symbols.radiation_field.symbol
        av = reac.symbols.visual_extinction.symbol

        habing = self.symbols.habing_field_photon_number.symbol
        crphot = self.symbols.cosmic_ray_induced_photon_number.symbol
        opt_uvd = self.symbols.photon_desorption_option.symbol
        cov = self.symbols.coverage.symbol
        nMono = self.symbols.active_monolayer_number.symbol
        garea = self.symbols.grain_surface_area.symbol

        sym_phot = f"{radfield}*{habing}*exp(-{av}*3.02) + {crphot} * ({crrate}/{zism})"

        spec = reac.reactants[0]
        rate = f"{opt_uvd} * {cov} * ({sym_phot}) * {spec.photon_yield()} * {nMono} * {garea}"
        return rate

    def rate_cosmicray_desorption(self, reac: Reaction) -> str:

        super().rate_cosmicray_desorption(reac)

        crrate = reac.symbols.cosmic_ray_ionization_rate.symbol
        zism = reac.symbols.ism_cosmic_ray_ionization_rate.symbol

        opt_crd = self.symbols.cosmic_ray_desorption_option.symbol
        cov = self.symbols.coverage.symbol
        nMono = self.symbols.active_monolayer_number.symbol
        sites = self.symbols.grain_surface_sites_density.symbol
        densites = self.symbols.surface_sites_density.symbol
        duty = self.symbols.cosmic_ray_desorption_duty_cycle.symbol
        Tcr = self.symbols.max_cosmic_ray_desorption_temperature.symbol

        spec = reac.reactants[0]
        rate = " * ".join(
            [
                f"{opt_crd} * {cov}",
                f"{duty} * {nMono} * {densites}",
                f"({crrate}/{zism})",
                f"sqrt(2.0*{sites}*kerg*eb_{spec.alias}/(pi*pi*amu*{spec.A}))",
                f"exp(-eb_{spec.alias}/{Tcr})",
            ]
        )
        return rate

    def rate_electron_capture(self, reac: Reaction) -> str:

        super().rate_electron_capture(reac)

        tgas = reac.symbols.temperature.symbol

        r = self.symbols.grain_radius.symbol

        rate = f"pi * {r} * {r} * sqrt(8.0*kerg*({tgas})/pi/amu/meu)"
        return rate

    def rate_recombination(self, reac: Reaction) -> str:

        super().rate_recombination(reac)

        tgas = reac.symbols.temperature.symbol

        r = self.symbols.grain_radius.symbol
        gdens = self.symbols.grain_density.symbol

        [spec] = [s for s in reac.reactants if not s.is_grain]
        a = reac.alpha

        rate = " * ".join(
            [
                f"{a} * pi * {r} * {r} * {gdens}",
                f"sqrt(8.0*kerg*{tgas}/(pi*amu*{spec.A}))",
                f"(1.0 + pow(echarge, 2.0)/{r}/kerg/{tgas})",
                f"(1.0 + sqrt(2.0*pow(echarge, 2.0)/({r}*kerg*{tgas}+2.0*pow(echarge, 2.0))))",
            ]
        )
        return rate

    def _rate_surface(self, reac: Reaction) -> str:

        re1, re2 = reac.reactants
        eb1, nmass1, eb2, nmass2 = re1.eb, re1.A, re2.eb, re2.A
        a = reac.alpha

        tdust = reac.symbols.dust_temperature.symbol

        gdens = self.symbols.grain_density.symbol
        freq = self.symbols.classical_diffusion_rate_factor.symbol
        quan = self.symbols.quantum_diffusion_rate_factor.symbol
        hop = self.symbols.surface_hopping_ratio.symbol
        cov = self.symbols.coverage.symbol
        nMono = self.symbols.active_monolayer_number.symbol
        unisites = self.symbols.unit_surface_sites.symbol
        densites = self.symbols.surface_sites_density.symbol

        afreq = f"{freq} * sqrt({eb1}/{nmass1})"
        adiff = f"{afreq} * exp(-{eb1}*{hop}/{tdust})/{unisites}"
        aquan = f"{afreq} * exp({quan} * sqrt({hop}*{nmass1}*{eb1})) / {unisites}"

        bfreq = f"{freq} * sqrt({eb2}/{nmass2})"
        bdiff = f"{bfreq} * exp(-{eb2}*{hop}/{tdust})/{unisites}"
        bquan = f"{bfreq} * exp({quan} * sqrt({hop}*{nmass2}*{eb2})) / {unisites}"

        kappa = f"exp(-{a}/{tdust})"
        kquan = f"exp({quan} * sqrt((({nmass1}*{nmass2})/({nmass1}+{nmass2}))*{a}))"

        rate = ""
        if re1.name in ["GH", "GH2"] and re2.name in ["GH", "GH2"]:
            rate = " * ".join(
                [
                    f"fmax({kappa}, {kquan})",
                    f"(fmax({adiff}, {aquan})+fmax({bdiff}, {bquan}))",
                    f"pow(({nMono}*{densites}), 2.0) / {gdens}",
                ]
            )
        elif re1.name in ["GH", "GH2"]:
            rate = " * ".join(
                [
                    f"fmax({kappa}, {kquan})",
                    f"(fmax({adiff}, {aquan})+{bdiff})",
                    f"pow(({nMono}*{densites}), 2.0) / {gdens}",
                ]
            )
        elif re2.name in ["GH", "GH2"]:
            rate = " * ".join(
                [
                    f"fmax({kappa}, {kquan})",
                    f"({adiff}+fmax({bdiff}, {bquan}))",
                    f"pow(({nMono}*{densites}), 2.0) / {gdens}",
                ]
            )
        else:
            rate = " * ".join(
                [
                    f"{kappa} * ({adiff}+{bdiff})",
                    f"pow(({nMono}*{densites}), 2.0) / {gdens}",
                ]
            )

        rate = " * ".join([rate, f"{cov}", f"{cov}"])
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
