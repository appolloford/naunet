import logging
from .. import settings
from ..species import Species
from .reaction import Reaction, ReactionType
from sympy.codegen.cfunctions import exp


class KIDAReaction(Reaction):
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
        zeta = settings.user_symbols["CRIR"]
        Tgas = settings.user_symbols["Temperature"]
        Av = settings.user_symbols["VisualExtinction"]
        if formula == 1:
            return a * zeta
        elif formula == 2:
            return a * exp(-c * Av)
        elif formula == 3:
            return a * (Tgas / 300.0) ** b * exp(-c / Tgas)
        else:
            raise RuntimeError(
                f"Formula {formula} has not been defined! Please extend the definition"
            )

    def _parse_string(self, react_string) -> None:
        react_string = react_string.strip()
        if react_string != "":
            rlen = 34  # length of the string containing reactants
            plen = 56  # length of the string containing products
            # print(react_string[:rlen].split())
            # print(react_string[rlen : rlen + plen].split())
            self.reactants = [
                self.create_species(r) for r in react_string[:rlen].split()
            ]
            self.products = [
                self.create_species(p) for p in react_string[rlen : rlen + plen].split()
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
            self.reaction_type = ReactionType(self.formula)
