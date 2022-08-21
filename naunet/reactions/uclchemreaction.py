from enum import IntEnum
from ..component import VariableType as vt
from ..grains.grain import Grain
from .reaction import Reaction
from .reactiontype import ReactionType as BasicType


class UCLCHEMReaction(Reaction):
    """
    The reaction format of UCLCHEM Makerates output.
    """

    class ReactionType(IntEnum):
        # two-body gas-phase reaction
        UCLCHEM_MA = BasicType.GAS_TWOBODY
        # direct cosmic-ray ionisation
        UCLCHEM_CR = BasicType.GAS_COSMICRAY
        # cosmic-ray-induced photoreaction
        UCLCHEM_CP = BasicType.GAS_UMIST_CRPHOT
        # photoreaction
        UCLCHEM_PH = BasicType.GAS_PHOTON
        # accretion
        UCLCHEM_FR = BasicType.GRAIN_FREEZE
        # thermal desorption
        UCLCHEM_TH = BasicType.GRAIN_DESORB_THERMAL
        # cosmic-ray-induced thermal desorption
        UCLCHEM_CD = BasicType.GRAIN_DESORB_COSMICRAY
        # photodesorption
        UCLCHEM_PD = BasicType.GRAIN_DESORB_PHOTON
        # reactive desorption
        UCLCHEM_RD = BasicType.GRAIN_DESORB_REACTIVE
        # H2-formation-induced desorption
        UCLCHEM_HD = BasicType.GRAIN_DESORB_H2
        # surface diffusion
        UCLCHEM_DF = BasicType.SURFACE_DIFFUSION

    reactant2type = {
        "CRP": ReactionType.UCLCHEM_CR,
        "PHOTON": ReactionType.UCLCHEM_PH,
        "CRPHOT": ReactionType.UCLCHEM_CP,
        "FREEZE": ReactionType.UCLCHEM_FR,
        "DESOH2": ReactionType.UCLCHEM_HD,
        "DESCR": ReactionType.UCLCHEM_CD,
        "DEUVCR": ReactionType.UCLCHEM_PD,
        "THERM": ReactionType.UCLCHEM_TH,
        "DIFF": ReactionType.UCLCHEM_DF,
        "CHEMDES": ReactionType.UCLCHEM_RD,
    }

    def __init__(self, react_string) -> None:
        super().__init__(format="uclchem", react_string=react_string)

        self.register("ism_cosmic_ray_ionization_rate", ("zism", 1.3e-17, vt.constant))
        self.register("x_ray_ionization_rate", ("zeta_xr", 0.0, vt.param))
        self.register("radiation_field", ("G0", 1.0, vt.param))
        self.register("dust_temperature", ("Tgas", None, vt.param))
        self.register("H2_column_density", ("h2col", "0.5*1.59e21*Av", vt.derived))
        self.register("CO_column_density", ("cocol", "1e-5 * h2col", vt.derived))
        self.register(
            "character_wavelength",
            ("lambdabar", "GetCharactWavelength(h2col, cocol)", vt.derived),
        )
        self.register(
            "H2_shielding_factor",
            (
                "H2shielding",
                "GetShieldingFactor(IDX_H2I, h2col, h2col, Tgas, 1)",
                vt.derived,
            ),
        )
        self.register(
            "H2_formation_rate",
            ("H2formation", "1.0e-17 * sqrt(Tgas) * nH", vt.derived),
        )
        self.register(
            "H2_dissociation_rate",
            (
                "H2dissociation",
                "5.1e-11 * G0 * GetGrainScattering(Av, 1000.0) * H2shielding",
                vt.derived,
            ),
        )

    def rateexpr(self, grain: Grain = None) -> str:
        a = self.alpha
        b = self.beta
        c = self.gamma
        rtype = self.reaction_type

        zeta = f"(zeta / zism)"  # uclchem has cosmic-ray ionization rate in unit of 1.3e-17s-1

        re1 = self.reactants[0]
        # re2 = self.reactants[1] if len(self.reactants) > 1 else None

        if grain:
            grain.sym_av = "Av"
            grain.sym_tgas = "Tgas"
            grain.sym_tdust = "Tgas"
            grain.sym_radfield = "G0"
            grain.sym_crrate = f"{zeta}"
            grain.sym_h2form = f"1.0e-17 * sqrt(Tgas) * y[IDX_HI] * nH"

        # two-body gas-phase reaction
        if rtype == self.ReactionType.UCLCHEM_MA:
            rate = " * ".join(
                s
                for s in [
                    f"{a}",
                    f"pow(Tgas/300.0, {b})" if b else "",
                    f"exp(-{c}/Tgas)" if c else "",
                ]
                if s
            )

        # direct cosmic-ray ionisation
        elif rtype == self.ReactionType.UCLCHEM_CR:
            rate = f"{a} * {zeta}"

        # cosmic-ray-induced photoreaction
        elif rtype == self.ReactionType.UCLCHEM_CP:
            rate = f"{a} * {zeta} * pow(Tgas/300.0, {b}) * {c} / (1.0 - omega)"

        # photoreaction
        elif rtype == self.ReactionType.UCLCHEM_PH:
            rate = f"G0 * {a} * exp(-{c}*Av) / 1.7"  # convert habing to Draine
            if re1.name in ["CO"]:
                shield = f"GetShieldingFactor(IDX_{re1.alias}, h2col, {re1.name.lower()}col, Tgas, 1)"
                rate = f"(2.0e-10) * G0 * {shield} * GetGrainScattering(Av, lambdabar) / 1.7"

        # accretion
        elif rtype == self.ReactionType.UCLCHEM_FR:
            rate = grain.rateexpr(self)

        # thermal desorption
        elif rtype == self.ReactionType.UCLCHEM_TH:
            rate = grain.rateexpr(self)

        # cosmic-ray-induced thermal desorption
        elif rtype == self.ReactionType.UCLCHEM_CD:
            rate = grain.rateexpr(self)

        # photodesorption
        elif rtype == self.ReactionType.UCLCHEM_PD:
            rate = grain.rateexpr(self)

        # H2 formation induced desorption
        elif rtype == self.ReactionType.UCLCHEM_HD:
            rate = grain.rateexpr(self)

        else:
            raise ValueError(f"Unsupported type: {rtype}")

        rate = self._beautify(rate)
        return rate

    def _parse_string(self, react_string) -> None:

        kwlist = [
            "NAN",
            "FREEZE",
            "DESOH2",
            "DESCR",
            "DEUVCR",
            "THERM",
            "DIFF",
            "CHEMDES",
        ]

        if react_string.strip() != "":

            *rpspec, a, b, c, lt, ut = react_string.split(",")

            self.reaction_type = self.reactant2type.get(
                rpspec[1], self.ReactionType.UCLCHEM_MA
            )

            # Turn off Freeze-out reaction beyond 30K
            if self.reaction_type == self.ReactionType.UCLCHEM_FR:
                lt, ut = 0, 30

            self.alpha = float(a)
            self.beta = float(b)
            self.gamma = float(c)
            self.temp_min = float(lt)
            self.temp_max = float(ut)
            # UCLCHEM program does not check the temperature range, the next two lines are used for benchmark test
            # self.temp_min = 10.0
            # self.temp_max = 41000.0

            reactants = [r for r in rpspec[0:3] if r not in kwlist]
            products = [p for p in rpspec[3:7] if p not in kwlist]

            self.reactants = [
                self._create_species(r) for r in reactants if self._create_species(r)
            ]
            self.products = [
                self._create_species(p) for p in products if self._create_species(p)
            ]
