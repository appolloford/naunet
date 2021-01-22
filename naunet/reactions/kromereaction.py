import logging
import re
from sympy import sympify
from sympy.codegen.cfunctions import exp
from ..species import Species
from ..settings import pseudo_element_list, user_symbols
from .reaction import Reaction, ReactionType


class KROMEReaction(Reaction):
    def __init__(self, react_string, format, *args, **kwargs) -> None:
        super().__init__(react_string)

        self.database = "krome"
        self.alpha = 0.0
        self.beta = 0.0
        self.gamma = 0.0
        self.kromeformat = format.lower().strip()
        self.rate_string = None

        self._parse_string(react_string)

    def rate_func(self):

        # print(self.rate_string)
        rate = re.sub(r"(\d\.?)d(\-?\d)", r"\1e\2", self.rate_string)
        return sympify(rate)

    def _parse_string(self, react_string, *argc, **kwargs) -> None:

        react_string = react_string.strip()
        if react_string != "" and react_string[0] != "#":
            kwords = self.kromeformat.split(",")

            for key, value in zip(kwords, react_string.split(",")):
                if value == "":
                    continue
                elif key == "r":
                    self.reactants.append(Species(value))
                elif key == "p":
                    self.products.append(Species(value))
                elif key == "tmin":
                    if value.upper() not in ["N", "NONE", "N/A", "NO", ""]:
                        for opstr in ["<", ">", ".LE.", ".GE.", ".LT.", ".GT."]:
                            value = value.replace(opstr, "")
                        value = value.replace("d", "e")
                        self.temp_min = float(value)
                elif key == "tmax":
                    if value.upper() not in ["N", "NONE", "N/A", "NO", ""]:
                        for opstr in ["<", ">", ".LE.", ".GE.", ".LT.", ".GT."]:
                            value = value.replace(opstr, "")
                        value = value.replace("d", "e")
                        self.temp_max = float(value)
                elif key == "rate":
                    self.rate_string = value.replace("dexp", "exp")
