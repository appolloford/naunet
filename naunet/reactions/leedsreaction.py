import logging
from .. import settings
from ..species import Species
from .reaction import Reaction, ReactionType


class LEEDSReaction(Reaction):
    """
    The reaction format used in Walsh et al. (2015) model. Named
    to LEEDSReaction because she is in University of Leeds now...
    The name can be changed anytime
    """

    def __init__(self, react_string, *args, **kwargs) -> None:
        super().__init__(react_string)

        self.database = "LEEDS"
        self.alpha = 0.0
        self.beta = 0.0
        self.gamma = 0.0
        self.rtype = None

        self._parse_string(react_string)

    def rate_func(self):
        a = self.alpha
        b = self.beta
        c = self.gamma
        rtype = self.rtype
        zeta = settings.user_symbols["CRIR"]
        Tgas = settings.user_symbols["Temperature"]
        Av = settings.user_symbols["VisualExtinction"]

        # TODO: finish the remaining type of reactions
        if rtype == 1:
            rate = f"{a} * pow({Tgas}/300.0, {b}) * exp(-{c}/{Tgas})"
        else:
            raise RuntimeError(
                f"Type {rtype} has not been defined! Please extend the definition"
            )

        rate = self._beautiy(rate)
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
                        Species(r)
                        for r in clip.split()
                        if r not in settings.pseudo_element_list
                    ]

                elif label == "prod":
                    self.products = [
                        Species(p)
                        for p in clip.split()
                        if p not in settings.pseudo_element_list
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
