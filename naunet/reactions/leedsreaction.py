import logging
from .. import settings
from ..dusts.dust import Dust
from .reaction import Reaction, ReactionType as BasicType


class LEEDSReaction(Reaction):
    """
    The reaction format used in Walsh et al. (2015) model. Named
    to LEEDSReaction because she is in University of Leeds now...
    The name can be changed anytime
    """

    consts = {
        "zism": 1.6e-3,
        "habing": 1e8,
        "crphot": 1e4,
    }
    vars = {
        "Hnuclei": "nH",
        "CRIR": "zeta_cr",
        "XRAY": "zeta_xr",
        "Temperature": "Tgas",
        "DustTemperature": "Tdust",
        "VisualExtinction": "Av",
        "UVPHOT": "uv",
        "G0": "G0",
    }
    user_var = []

    def __init__(self, react_string, *args, dust: Dust = None, **kwargs) -> None:
        super().__init__(react_string)

        self.database = "LEEDS"
        self.alpha = 0.0
        self.beta = 0.0
        self.gamma = 0.0
        self.rtype = None
        self.dust = dust

        settings.surface_symbol = "G"

        self._parse_string(react_string)

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
            rg = self.dust.vars.get("Radius")
            albedo = self.dust.vars.get("Albedo")
            sites = self.dust.vars.get("SurfaceSites")
            nmono = self.dust.vars.get("MonoLayers")
            duty = self.dust.vars.get("DutyCycle")
            Tcr = self.dust.vars.get("CRDesorptionTemperature")
            gdens = "(y[IDX_GRAIN0I]+y[IDX_GRAINM])"
            garea = f"(4*pi*{rg}*{rg}) * {gdens}"
            densites = f"{sites} * (4*pi*{rg}*{rg}) * {gdens}"

        # TODO: finish the remaining type of reactions
        if rtype == 1:
            rate = f"{a} * pow({Tgas}/300.0, {b}) * exp(-{c}/{Tgas})"
        elif rtype == 2:
            rate = f"{a} * ({cr} + {xr}) / {zism}"
        elif rtype == 3:
            rate = f"{a} * pow({Tgas}/300.0, {b}) * {c} * ({cr} + {xr}) / {zism} / (1.0 - {albedo})"
        elif rtype == 4:
            # TODO: shielding
            rate = f"{G0} * {a} * exp(-{c}*{Av})"
        elif rtype == 5:
            # TODO: X-ray
            rate = "0.0"
        elif rtype == 6:
            rate = f"{a} * pi * pow({rg}, 2) * {gdens} * sqrt(8.0*kerg*{Tgas}/pi/amu/{c}) * (1.0+pow(echarge, 2)/{rg}/kerg/Tgas) * (1.0 + sqrt(2.0*pow(echarge, 2)/({rg}*kerg*{Tgas}+2.0*pow(echarge, 2))))"
        elif rtype == 7:
            rate = f"{a} * pi * pow({rg}, 2) * {gdens} * sqrt(8.0 * kerg * {Tgas}/ (pi*amu*{c}))"
        elif rtype == 8:
            rate = f"sqrt(2.0*{sites}*kerg*eb_{re1.alias}/(pi*pi*amu*{c})) * {nmono} * {densites} * exp(-eb_{re1.alias}/{Tdust})"
        elif rtype == 9:
            rate = f"({cr}/{zism}) * {duty} * sqrt(2.0*{sites}*kerg*eb_{re1.alias}/(pi*pi*amu*{c})) * {nmono} * {densites} * exp(-eb_{re1.alias}/{Tcr})"
        elif rtype == 10:
            uvphot = f"{G0}*{habing}*exp(-{Av}*3.02) + {crphot} * {cr}/{zism}"
            rate = f"{uvphot} * {re1.photon_yield} * {nmono} * {garea}"
        elif rtype == 11:
            rate = f"{a} * ({xr}+{cr})/{zism} * pow({Tgas}/300.0, {b}) * {c} / (1.0 - {albedo})"
        elif rtype == 12:
            # TODO: shielding
            rate = f"{G0} * {a} * exp(-{c}*{Av})"
        elif rtype == 13:
            # TODO: diffusion
            rate = "0.0"
        elif rtype == 14:
            # TODO: reaction desorption
            rate = "0.0"
        elif rtype in range(15, 20):
            rate = "0.0"
        elif rtype == 20:
            rate = f"pi * {rg} * {rg} * sqrt(8.0*kerg*{Tgas}/pi/amu/me)"
        else:
            raise RuntimeError(
                f"Type {rtype} has not been defined! Please extend the definition"
            )

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
                else:
                    raise RuntimeError("Unknown label inserted! Please check")

                stidx += len
