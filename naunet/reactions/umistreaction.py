from enum import IntEnum
from ..grains.grain import Grain
from .reaction import Reaction
from ..reactiontype import ReactionType as BasicType


class UMISTReaction(Reaction):
    format = "umist"

    class ReactionType(IntEnum):
        # A representation for other sub categories
        UMIST_TWOBODY = BasicType.GAS_TWOBODY
        UMIST_AD = BasicType.GAS_TWOBODY  # Associative Detachment
        UMIST_CD = BasicType.GAS_TWOBODY  # Collisional Dissociation
        UMIST_CE = BasicType.GAS_TWOBODY  # Charge Exchange
        UMIST_CP = BasicType.GAS_COSMICRAY  # Cosmic-Ray Proton (CRP)
        UMIST_CR = BasicType.GAS_UMIST_CRPHOT  # Cosmic-Ray Photon (CRPHOT)
        UMIST_DR = BasicType.GAS_TWOBODY  # Dissociative Recombination
        UMIST_IN = BasicType.GAS_TWOBODY  # Ion-Nuetral
        UMIST_MN = BasicType.GAS_TWOBODY  # Mutual Neutralisation
        UMIST_NN = BasicType.GAS_TWOBODY  # Nuetral-Neutral
        UMIST_PH = BasicType.GAS_PHOTON  # Photoprocess
        UMIST_RA = BasicType.GAS_TWOBODY  # Radiative Association
        UMIST_REA = BasicType.GAS_TWOBODY  # Radiative Electron Attachment
        UMIST_RR = BasicType.GAS_TWOBODY  # Radiative Recombination

    # map the UMIST code reaction types, see McElroy+2013
    code2type = {
        "AD": ReactionType.UMIST_AD,
        "CD": ReactionType.UMIST_CD,
        "CE": ReactionType.UMIST_CE,
        "CP": ReactionType.UMIST_CP,
        "CR": ReactionType.UMIST_CR,
        "DR": ReactionType.UMIST_DR,
        "IN": ReactionType.UMIST_IN,
        "MN": ReactionType.UMIST_MN,
        "NN": ReactionType.UMIST_NN,
        "PH": ReactionType.UMIST_PH,
        "RA": ReactionType.UMIST_RA,
        "REA": ReactionType.UMIST_REA,
        "RR": ReactionType.UMIST_RR,
    }

    def __init__(self, react_string: str) -> None:
        # Extra attributes in UMISTReaction
        self.code = None

        super().__init__(react_string=react_string)

    def rateexpr(self, grain: Grain = None) -> str:
        a = self.alpha
        b = self.beta
        c = self.gamma
        rtype = self.reaction_type

        if rtype == self.ReactionType.UMIST_TWOBODY:
            rate = " * ".join(
                s
                for s in [
                    f"{a}",
                    f"pow(Tgas/300.0, {b})" if b else "",
                    f"exp(-{c}/Tgas)" if c else "",
                ]
                if s
            )
        elif rtype == self.ReactionType.UMIST_PH:
            rate = f"{a} * exp(-{c}*Av)"
        elif rtype == self.ReactionType.UMIST_CP:
            rate = f"{a}"
        elif rtype == self.ReactionType.UMIST_CR:
            rate = f"{a} * pow(Tgas/300.0, {b}) * {c} / (1-omega)"
        else:
            raise RuntimeError(
                f"Code {self.code} has not been defined! Please extend the definition"
            )

        rate = self._beautify(rate)
        return rate

    def _parse_string(self, react_string) -> None:
        self.source = "umist"

        react_string = react_string.strip()
        if react_string != "":
            idx, code, *rps, _, a, b, c, lt, ut = react_string.split(":")[:14]
            # print(id, rps)
            self.reactants = [
                self._create_species(r) for r in rps[0:2] if self._create_species(r)
            ]
            self.products = [
                self._create_species(p) for p in rps[2:6] if self._create_species(p)
            ]

            self.alpha = float(a)
            self.beta = float(b)
            self.gamma = float(c)
            self.temp_min = float(lt)
            self.temp_max = float(ut)
            self.idxfromfile = int(idx)
            self.code = code
            self.reaction_type = self.code2type.get(self.code)
