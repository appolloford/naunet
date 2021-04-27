import logging
from enum import IntEnum
from .reaction import Reaction, ReactionType as BasicType


class KIDAReaction(Reaction):
    class ReactionType(IntEnum):
        KIDA_MA = BasicType.GAS_TWOBODY  # Modified Arrhenius
        KIDA_CR = BasicType.GAS_COSMICRAY  # Cosmic-ray ionization
        KIDA_PD = BasicType.GAS_PHOTON  # Photo-dissociation (Draine)
        KIDA_TB = BasicType.GAS_THREEBODY  # Three-body
        KIDA_IP1 = BasicType.GAS_KIDA_IP1  # ionpol1
        KIDA_IP2 = BasicType.GAS_KIDA_IP2  # ionpol2

    # map the formula to KIDA reaction types
    formula2type = {
        1: ReactionType.KIDA_CR,
        2: ReactionType.KIDA_PD,
        3: ReactionType.KIDA_MA,
        4: ReactionType.KIDA_IP1,
        5: ReactionType.KIDA_IP2,
        6: ReactionType.KIDA_TB,
    }

    vars = {
        "Hnuclei": "nH",
        "CRIR": "zeta",
        "Temperature": "Tgas",
        "VisualExtinction": "Av",
        "UVPHOT": "uv",
    }
    user_var = []

    def __init__(self, react_string, *args, **kwargs) -> None:
        super().__init__(react_string)

        self.database = "KIDA"
        self.alpha = 0.0
        self.beta = 0.0
        self.gamma = 0.0
        self.formula = -1
        self.itype = -1

        self._parse_string(react_string)

    def rate_func(self):
        a = self.alpha
        b = self.beta
        c = self.gamma
        formula = self.formula
        zeta = KIDAReaction.vars["CRIR"]
        Tgas = KIDAReaction.vars["Temperature"]
        Av = KIDAReaction.vars["VisualExtinction"]
        if formula == 1:
            rate = f"{a} * {zeta}"
        elif formula == 2:
            rate = f"{a} * exp(-{c}*{Av})"
        elif formula == 3:
            rate = f"{a} * pow({Tgas}/300.0, {b}) * exp(-{c}/{Tgas})"
        elif formula == 4:
            rate = f"{a} * {b} * (0.62 + 0.4767*{c}*sqrt(300.0/{Tgas}))"
        elif formula == 5:
            rate = f"{a} * {b} * (1 + 0.0967*{c}*sqrt(300.0/{Tgas}) + c*c*(300.0/{Tgas})/10.526)"
        elif formula == 6:
            raise NotImplementedError("Three-body reactions formula is not implemented")
        else:
            raise RuntimeError(
                f"Formula {formula} has not been defined! Please extend the definition"
            )

        rate = self._beautify(rate)
        return rate

    def _parse_string(self, react_string) -> None:
        react_string = react_string.strip()
        if react_string != "":
            rlen = 34  # length of the string containing reactants
            plen = 56  # length of the string containing products
            # print(react_string[:rlen].split())
            # print(react_string[rlen : rlen + plen].split())
            self.reactants = [
                self.create_species(r)
                for r in react_string[:rlen].split()
                if self.create_species(r)
            ]
            self.products = [
                self.create_species(p)
                for p in react_string[rlen : rlen + plen].split()
                if self.create_species(p)
            ]

            a, b, c, _, _, _, itype, lt, ut, form, _, _, _ = react_string[
                rlen + plen :
            ].split()

            self.alpha = float(a)
            self.beta = float(b)
            self.gamma = float(c)
            self.itype = int(itype)
            self.temp_min = float(lt)
            self.temp_max = float(ut)
            self.formula = int(form)
            if self.formula < 1 or self.formula > 6:
                logging.warning(
                    f"Formula {form} is not valid in reaction {self}, change to formula = 3."
                )
                self.formula = 3
            self.reaction_type = self.formula2type.get(self.formula)
