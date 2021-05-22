import logging
from enum import IntEnum
from ..species import Species
from ..dusts.dust import Dust
from .reaction import Reaction, ReactionType as BasicType


class LEEDSReaction(Reaction):
    """
    The reaction format used in Walsh et al. (2015) model. Named
    to LEEDSReaction because she is in University of Leeds now...
    The name can be changed anytime
    """

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
        LEEDS_XR = BasicType.GAS_LEEDS_XRAY
        # cation-grain recombination
        LEEDS_RC = BasicType.GAS_LEEDS_RECOM
        # accretion
        LEEDS_FR = BasicType.GRAIN_FREEZE
        # thermal desorption
        LEEDS_TH = BasicType.GRAIN_DESORPT_THERMAL
        # cosmic-ray-induced thermal desorption
        LEEDS_CD = BasicType.GRAIN_DESORPT_COSMICRAY
        # photodesorption
        LEEDS_PD = BasicType.GRAIN_DESORPT_PHOTON
        # grain-surface cosmic-ray-induced photoreaction
        LEEDS_SC = BasicType.SURFACE_COSMICRAY
        # grain-surface photoreaction
        LEEDS_SP = BasicType.SURFACE_PHOTON
        # two-body grain-surface reaction
        LEEDS_SB = BasicType.SURFACE_TWOBODY
        # reactive desorption
        LEEDS_RD = BasicType.GRAIN_DESORPT_REACTIVE
        # grain electron capture rate
        LEEDS_EC = BasicType.GAS_LEEDS_ECAPTURE

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

    consts = {
        "gism": 1.6e-3,
        "habing": 1e8,
        "zism": 1.3e-17,
        "crphot": 1e4,
        "hbar": 1.054571726e-27,
    }
    vars = {
        "Hnuclei": "nH",
        "CRIR": "zeta_cr",
        "XRAY": "zeta_xr",
        "Temperature": "Tgas",
        "DustTemperature": "Tdust",
        "VisualExtinction": "Av",
        "G0": "G0",
    }
    user_var = [
        "double h2col = 0.5*1.59e21*Av",
        "double cocol = 1e-5 * h2col",
        "double n2col = 1e-5 * h2col",
    ]

    def __init__(self, react_string, *args, dust: Dust = None, **kwargs) -> None:
        super().__init__(react_string)

        self.database = "LEEDS"
        self.alpha = 0.0
        self.beta = 0.0
        self.gamma = 0.0
        self.rtype = None
        self.dust = dust

        self._parse_string(react_string)

    @classmethod
    def initialize(cls) -> None:
        Species.surface_prefix = "G"

    @classmethod
    def finalize(cls) -> None:
        Species.surface_prefix = "#"

    def rate_func(self):
        a = self.alpha
        b = self.beta
        c = self.gamma
        rtype = self.rtype

        cr = LEEDSReaction.vars.get("CRIR")
        xr = LEEDSReaction.vars.get("XRAY")
        Tgas = LEEDSReaction.vars.get("Temperature")
        Tdust = LEEDSReaction.vars.get("DustTemperature")
        Av = LEEDSReaction.vars.get("VisualExtinction")
        G0 = LEEDSReaction.vars.get("G0")
        zism = LEEDSReaction.consts.get("zism")
        habing = LEEDSReaction.consts.get("habing")
        crphot = LEEDSReaction.consts.get("crphot")

        if self.dust:
            re1 = self.reactants[0]
            if self.rtype in [13, 14]:
                re2 = self.reactants[1]
            rg = self.dust.vars.get("Radius")
            albedo = self.dust.vars.get("Albedo")
            sites = self.dust.vars.get("SurfaceSites")
            hop = self.dust.vars.get("HOPRatio")
            nmono = self.dust.vars.get("MonoLayers")
            duty = self.dust.vars.get("DutyCycle")
            Tcr = self.dust.vars.get("CRDesorptionTemperature")
            branch = self.dust.vars.get("BranchRatio")

        if rtype == 1:
            rate = f"{a} * pow({Tgas}/300.0, {b}) * exp(-{c}/{Tgas})"
        elif rtype == 2:
            rate = f"{a} * ({cr} + {xr}) / {zism}"
        elif rtype == 3:
            rate = f"{a} * (({cr} + {xr}) / {zism}) * pow({Tgas}/300.0, {b}) * {c} / (1.0 - {albedo})"
        elif rtype == 4:
            rate = f"{G0} * {a} * exp(-{c}*{Av})"
            if re1.name in ["H2", "CO", "N2"]:
                shield = f"shieldingfactor(IDX_{re1.alias}, h2col, {re1.name.lower()}col, {Tgas}, 0)"
                rate = f"{rate} * {shield}"
        elif rtype == 5:
            # TODO: X-ray
            rate = "0.0"
        elif rtype == 6:
            rate = f"{a} * pi * pow({rg}, 2.0) * gdens * sqrt(8.0*kerg*{Tgas}/(pi*amu*{c})) * (1.0+pow(echarge, 2)/{rg}/kerg/{Tgas}) * (1.0 + sqrt(2.0*pow(echarge, 2.0)/({rg}*kerg*{Tgas}+2.0*pow(echarge, 2.0))))"
        elif rtype == 7:
            rate = f"{a} * pi * pow({rg}, 2.0) * gdens * sqrt(8.0 * kerg * {Tgas}/ (pi*amu*{c}))"
        elif rtype == 8:
            rate = f"sqrt(2.0*{sites}*kerg*eb_{re1.alias}/(pi*pi*amu*{c})) * {nmono} * densites * exp(-eb_{re1.alias}/{Tdust})"
        elif rtype == 9:
            rate = f"({cr}/{zism}) * {duty} * sqrt(2.0*{sites}*kerg*eb_{re1.alias}/(pi*pi*amu*{c})) * {nmono} * densites * exp(-eb_{re1.alias}/{Tcr})"
            rate = "0.0"
        elif rtype == 10:
            uvphot = f"{G0}*{habing}*exp(-{Av}*3.02) + {crphot} * {cr}/{zism}"
            rate = f"({uvphot}) * {re1.photon_yield} * {nmono} * garea"
        elif rtype == 11:
            rate = f"{a} * ({xr}+{cr})/{zism} * pow({Tgas}/300.0, {b}) * {c} / (1.0 - {albedo})"
        elif rtype == 12:
            rate = f"{G0} * {a} * exp(-{c}*{Av})"
            if re1.name in ["GH2", "GCO", "GN2"]:
                shield = f"shieldingfactor(IDX_{re1.alias[1:]}, h2col, {re1.name[1:].lower()}col, {Tgas}, 0)"
                rate = f"{rate} * {shield}"
        elif rtype == 13:
            afreq = f"freq * sqrt({re1.binding_energy}/{re1.massnumber})"
            adiff = f"{afreq} * exp(-{re1.binding_energy}*{hop}/{Tdust})/unisites"
            aquan = f"{afreq} * exp(quan * sqrt({hop}*{re1.massnumber}*{re1.binding_energy})) / unisites"

            bfreq = f"freq * sqrt({re2.binding_energy}/{re2.massnumber})"
            bdiff = f"{bfreq} * exp(-{re2.binding_energy}*{hop}/{Tdust})/unisites"
            bquan = f"{bfreq} * exp(quan * sqrt({hop}*{re2.massnumber}*{re2.binding_energy})) / unisites"

            kappa = f"exp(-{a}/{Tdust})"
            kquan = f"exp(quan * sqrt((({re1.massnumber}*{re2.massnumber})/({re1.massnumber}+{re2.massnumber}))*{a}))"

            if re1.name in ["GH", "GH2"] and re2.name in ["GH", "GH2"]:
                rate = f"fmax({kappa}, {kquan})*(fmax({adiff}, {aquan})+fmax({bdiff}, {bquan}))*pow(({nmono}*densites), 2.0)/gdens"
            elif re1.name in ["GH", "GH2"]:
                rate = f"fmax({kappa}, {kquan})*(fmax({adiff}, {aquan})+{bdiff})*pow(({nmono}*densites), 2.0)/gdens"
            elif re2.name in ["GH", "GH2"]:
                rate = f"fmax({kappa}, {kquan})*({adiff}+fmax({bdiff}, {bquan}))*pow(({nmono}*densites), 2.0)/gdens"
            else:
                rate = f"{kappa}*({adiff}+{bdiff})*pow(({nmono}*densites), 2.0)/gdens"
        elif rtype == 14:

            afreq = f"freq * sqrt({re1.binding_energy}/{re1.massnumber})"
            adiff = f"{afreq} * exp(-{re1.binding_energy}*{hop}/{Tdust})/unisites"
            aquan = f"{afreq} * exp(quan * sqrt({hop}*{re1.massnumber}*{re1.binding_energy})) / unisites"

            bfreq = f"freq * sqrt({re2.binding_energy}/{re2.massnumber})"
            bdiff = f"{bfreq} * exp(-{re2.binding_energy}*{hop}/{Tdust})/unisites"
            bquan = f"{bfreq} * exp(quan * sqrt({hop}*{re2.massnumber}*{re2.binding_energy})) / unisites"

            kappa = f"exp(-{a}/{Tdust})"
            kquan = f"exp(quan * sqrt((({re1.massnumber}*{re2.massnumber})/({re1.massnumber}+{re2.massnumber}))*{a}))"

            if re1.name in ["GH", "GH2"] and re2.name in ["GH", "GH2"]:
                rate = f"fmax({kappa}, {kquan})*(fmax({adiff}, {aquan})+fmax({bdiff}, {bquan}))*pow(({nmono}*densites), 2.0)/gdens"
            elif re1.name in ["GH", "GH2"]:
                rate = f"fmax({kappa}, {kquan})*(fmax({adiff}, {aquan})+{bdiff})*pow(({nmono}*densites), 2.0)/gdens"
            elif re2.name in ["GH", "GH2"]:
                rate = f"fmax({kappa}, {kquan})*({adiff}+fmax({bdiff}, {bquan}))*pow(({nmono}*densites), 2.0)/gdens"
            else:
                rate = f"{kappa}*({adiff}+{bdiff})*pow(({nmono}*densites), 2.0)/gdens"
            rate = f"{branch} * {rate}"
        elif rtype in range(15, 20):
            rate = "0.0"
        elif rtype == 20:
            rate = f"pi * {rg} * {rg} * sqrt(8.0*kerg*{Tgas}/pi/amu/meu)"
        else:
            raise RuntimeError(
                f"Type {rtype} has not been defined! Please extend the definition"
            )

        nsurface = len([r for r in self.reactants if r.is_surface])
        rate = rate + "*cov" * nsurface if rtype not in [11, 12] else rate

        rate = self._beautify(rate)
        return rate

    def _parse_string(self, react_string) -> None:
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
                    pass
                elif label == "reac":
                    self.reactants = [
                        self.create_species(r.replace("YC", "CH2OHC"))
                        for r in clip.split()
                        if self.create_species(r.replace("YC", "CH2OHC"))
                    ]

                elif label == "prod":
                    self.products = [
                        self.create_species(p.replace("YC", "CH2OHC"))
                        for p in clip.split()
                        if self.create_species(p.replace("YC", "CH2OHC"))
                    ]

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
