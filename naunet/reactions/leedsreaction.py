import logging
from enum import IntEnum
from ..component import VariableType as vt
from ..species import Species
from ..grains.grain import Grain
from .reaction import Reaction
from ..reactiontype import ReactionType as BasicType


class LEEDSReaction(Reaction):
    """
    The reaction format used in Walsh et al. (2015) model. Named
    to LEEDSReaction because she is in University of Leeds now...
    The name can be changed anytime
    """

    format = "leeds"

    class ReactionType(IntEnum):
        # two-body gas-phase reaction
        LEEDS_MA = BasicType.GAS_TWOBODY
        # direct cosmic-ray ionisation
        LEEDS_CR = BasicType.GAS_COSMICRAY
        # cosmic-ray-induced photoreaction
        LEEDS_CP = BasicType.GAS_UMIST_CRPHOT
        # photoreaction
        LEEDS_PH = BasicType.GAS_PHOTON
        # direct X-ray ionisation
        LEEDS_XR = BasicType.GAS_XRAY
        # cation-grain recombination
        LEEDS_RC = BasicType.GRAIN_RECOMINE
        # accretion
        LEEDS_FR = BasicType.GRAIN_FREEZE
        # thermal desorption
        LEEDS_TH = BasicType.GRAIN_DESORB_THERMAL
        # cosmic-ray-induced thermal desorption
        LEEDS_CD = BasicType.GRAIN_DESORB_COSMICRAY
        # photodesorption
        LEEDS_PD = BasicType.GRAIN_DESORB_PHOTON
        # grain-surface cosmic-ray-induced photoreaction
        LEEDS_SC = BasicType.SURFACE_COSMICRAY
        # grain-surface photoreaction
        LEEDS_SP = BasicType.SURFACE_PHOTON
        # two-body grain-surface reaction
        LEEDS_SB = BasicType.SURFACE_TWOBODY
        # reactive desorption
        LEEDS_RD = BasicType.GRAIN_DESORB_REACTIVE
        # grain electron capture rate
        LEEDS_EC = BasicType.GRAIN_ECAPTURE

    rtype2type = {
        1: ReactionType.LEEDS_MA,
        2: ReactionType.LEEDS_CR,
        3: ReactionType.LEEDS_CP,
        4: ReactionType.LEEDS_PH,
        5: ReactionType.LEEDS_XR,
        6: ReactionType.LEEDS_RC,
        7: ReactionType.LEEDS_FR,
        8: ReactionType.LEEDS_TH,
        9: ReactionType.LEEDS_CD,
        10: ReactionType.LEEDS_PD,
        11: ReactionType.LEEDS_SC,
        12: ReactionType.LEEDS_SP,
        13: ReactionType.LEEDS_SB,
        14: ReactionType.LEEDS_RD,
        20: ReactionType.LEEDS_EC,
    }

    def __init__(self, react_string: str) -> None:
        # Extra attributes in LEEDSReaction
        self.rtype = None

        super().__init__(react_string=react_string)

        self.register("ism_radiation_field", ("gism", 1.6e-3, vt.constant))
        self.register("ism_cosmic_ray_ionization_rate", ("zism", 1.3e-17, vt.constant))
        self.register(
            "cosmic_ray_ionization_rate",
            ("zeta_cr", 1.3e-17, vt.param),
            force_overwrite=True,
        )
        self.register("x_ray_ionization_rate", ("zeta_xr", 0.0, vt.param))
        self.register("radiation_field", ("G0", 1.0, vt.param))
        self.register(
            "dust_temperature", ("Tdust", None, vt.param), force_overwrite=True
        )
        self.register("H2_column_density", ("h2col", "0.5*1.59e21*Av", vt.derived))
        self.register("CO_column_density", ("cocol", "1e-5 * h2col", vt.derived))
        self.register("N2_column_density", ("n2col", "1e-5 * h2col", vt.derived))

        self.register(
            "sticking_coefficient1",
            (
                f"stick1",
                f"(1.0 / (1.0 + 4.2e-2*sqrt(Tgas+Tdust) + 2.3e-3*Tgas - 1.3e-7*Tgas*Tgas))",
                vt.derived,
            ),
        )
        self.register(
            "sticking_coefficient2",
            (
                f"stick2",
                f"exp(-1741.0/Tgas) / (1.0 + 5e-2*sqrt(Tgas+Tdust) + 1e-14*pow(Tgas, 4.0))",
                vt.derived,
            ),
        )
        self.register(
            "total_sticking_coefficient",
            (
                f"stick",
                f"stick1 + stick2",
                vt.derived,
            ),
        )

    def rateexpr(self, grain: Grain = None) -> str:
        a = self.alpha
        b = self.beta
        c = self.gamma
        rtype = self.rtype

        re1 = self.reactants[0]
        re2 = self.reactants[1] if len(self.reactants) > 1 else None

        # two-body gas-phase reaction
        if rtype == 1:
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
        elif rtype == 2:
            rate = f"{a} * (zeta_cr + zeta_xr) / zism"

        # cosmic-ray-induced photoreaction
        elif rtype == 3:
            rate = f"{a} * ((zeta_cr + zeta_xr) / zism) * pow(Tgas/300.0, {b}) * {c} / (1.0 - omega)"

        # photoreaction
        elif rtype == 4:
            rate = f"G0 * {a} * exp(-{c}*Av)"
            if re1.name in ["H2", "CO", "N2"]:
                shield = f"GetShieldingFactor(IDX_{re1.alias}, h2col, {re1.name.lower()}col, Tgas, 0)"
                rate = f"{rate} * {shield}"

        # TODO:
        # direct X-ray ionisation
        elif rtype == 5:
            rate = "0.0"

        # cation-grain recombination
        elif rtype == 6:
            rate = grain.rateexpr(self)

        # accretion
        elif rtype == 7:
            rate = grain.rateexpr(self)

        # thermal desorption
        elif rtype == 8:
            rate = grain.rateexpr(self)

        # cosmic-ray-induced thermal desorption
        elif rtype == 9:
            rate = grain.rateexpr(self)

        # photodesorption
        elif rtype == 10:
            rate = grain.rateexpr(self)

        # grain-surface cosmic-ray-induced photoreaction
        elif rtype == 11:
            zeta = f"(zeta_xr+zeta_cr)/zism"
            rate = f"{a} * ({zeta}) * pow(Tgas/300.0, {b}) * {c} / (1.0 - omega)"

        # grain-surface photoreaction
        elif rtype == 12:
            rate = f"G0 * {a} * exp(-{c}*Av)"
            if re1.name in ["GH2", "GCO", "GN2"]:
                spidx = f"IDX_{re1.alias[1:]}"
                coldens = f"{re1.name[1:].lower()}col"
                shield = f"GetShieldingFactor({spidx}, h2col, {coldens}, Tgas, 0)"
                rate = f"{rate} * {shield}"

        # two-body grain-surface reaction
        elif rtype == 13:
            rate = grain.rateexpr(self)

        # reactive desorption
        elif rtype == 14:
            rate = grain.rateexpr(self)

        # three-body association *
        # collisional dissociation *
        # collisional de-excitation of H2* *
        # Lyman-alpha photoreaction *
        # radiative de-excitation of H2* *
        elif rtype in range(15, 20):
            rate = "0.0"

        # grain electron capture rate
        elif rtype == 20:
            rate = grain.rateexpr(self)

        else:
            raise RuntimeError(
                f"Type {rtype} has not been defined! Please extend the definition"
            )

        rate = self._beautify(rate)
        return rate

    def _parse_string(self, react_string) -> None:
        self.source = "leeds"

        list_label = [
            "idx",
            "reac",
            "prod",
            "a",
            "b",
            "c",
            "lt",
            "ht",
            "type",
        ]
        list_strlen = [
            5,
            30,
            50,
            8,
            9,
            10,
            5,
            5,
            3,
        ]
        # react_string = react_string.strip()
        if react_string.strip() != "":
            stidx = 0
            for label, len in zip(list_label, list_strlen):
                clip = react_string[stidx : stidx + len]
                if label == "idx":
                    self.idxfromfile = int(clip)
                elif label == "reac":
                    reactants = [
                        self._create_species(
                            r.replace("YC", "CH2OHC"), surface_prefix="G"
                        )
                        for r in clip.split()
                    ]
                    # remove `None` values in the list
                    self.reactants = [s for s in reactants if s]

                elif label == "prod":
                    products = [
                        self._create_species(
                            p.replace("YC", "CH2OHC"), surface_prefix="G"
                        )
                        for p in clip.split()
                    ]
                    # remove `None` values in the list
                    self.products = [s for s in products if s]

                elif label == "a":
                    self.alpha = float(clip)
                elif label == "b":
                    self.beta = float(clip)
                elif label == "c":
                    self.gamma = float(clip)
                elif label == "lt":
                    self.temp_min = float(clip)
                elif label == "ht":
                    self.temp_max = float(clip)
                elif label == "type":
                    self.rtype = int(clip[1:])  # the first char is not used
                    self.reaction_type = self.rtype2type.get(self.rtype)
                else:
                    raise RuntimeError("Unknown label inserted! Please check")

                stidx += len
