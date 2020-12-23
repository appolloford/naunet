from gac.species import Species
from gac.reaction import Reaction, ReactionType
from gac.singleton import pseudo_element_list, user_symbols
from sympy.codegen.cfunctions import exp


class KIDAReaction(Reaction):
    def __init__(self, react_string) -> None:
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
        zeta = user_symbols["CRIR"]
        Tgas = user_symbols["Temperature"]
        Av = user_symbols["VisualExtinction"]
        if self.reaction_type == ReactionType.KIDA_CR:
            return a * zeta
        elif self.reaction_type == ReactionType.KIDA_PD:
            return a * exp(-c * Av)
        elif self.reaction_type == ReactionType.KIDA_MA:
            return a * (Tgas / 300.0) ** b * exp(-c / Tgas)
        else:
            raise RuntimeError(
                "Reaction rate formula could not be generated for reaction type: ".format(
                    ReactionType(self.reaction_type)
                )
            )

    def _parse_string(self, react_string) -> None:
        react_string = react_string.strip()
        if react_string != "":
            rlen = 34  # length of the string containing reactants
            plen = 56  # length of the string containing products
            # print(react_string[:rlen].split())
            # print(react_string[rlen : rlen + plen].split())
            self.reactants = [
                Species(r)
                for r in react_string[:rlen].split()
                if r not in pseudo_element_list
            ]
            self.products = [
                Species(p)
                for p in react_string[rlen : rlen + plen].split()
                if p not in pseudo_element_list
            ]
            a, b, c, _, _, _, itype, lt, ut, form, _, _, _ = react_string[
                rlen + plen :
            ].split()

            print(form)
            self.alpha = float(a)
            self.beta = float(b)
            self.gamma = float(c)
            self.itype = int(itype)
            self.temp_min = float(lt)
            self.temp_max = float(ut)
            self.formula = int(form)
            self.reaction_type = ReactionType(self.formula)
