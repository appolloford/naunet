import re
from ..grains.grain import Grain
from ..species import Species
from .reaction import Reaction
from .converter import ExpressionConverter


class KROMEReaction(Reaction):

    format = "krome"

    _kromerateconverter = ExpressionConverter("Fortran")

    def __init__(self, react_string: str) -> None:
        super().__init__(react_string=react_string)

    @classmethod
    def initialize(cls) -> None:
        cls.reacformat = "idx,r,r,r,p,p,p,p,tmin,tmax,rate"
        cls.varis = {
            "nH": None,
            "Tgas": None,
        }
        cls.locvars = []

    @classmethod
    def preprocessing(cls, line: str) -> str:
        if line.startswith(("#", "//")):
            return ""
        elif line.startswith("@format:"):
            cls.reacformat = line.replace("@format:", "")
            return ""
        elif line.startswith("@var"):
            if "Hnuclei" not in line:
                cls.locvars.append(line.replace("@var:", "").strip())
            return ""
        elif line.startswith("@common:"):
            commonlist = line.replace("@common:", "").strip().split(",")
            cls.varis.update(zip(commonlist, [None] * len(commonlist)))
            return ""
        else:
            return line.strip()

    def rateexpr(self, grain: Grain = None) -> str:

        rate = re.sub(r"(\d\.?)d(\-?\d)", r"\1e\2", self.rate_string)
        rate = re.sub(r"(idx_.?)p", r"\1II", rate)
        rate = re.sub(r"(idx_.?)m", r"\1M", rate)
        rate = re.sub(r"(idx_.?)\)", r"\1I)", rate)
        rate = rate.replace("Hnuclei", "nH")
        self._kromerateconverter.read(rate)
        rate = f"{self._kromerateconverter:c}"
        return rate

    def _parse_string(self, react_string) -> None:

        # Extra attributes in KROMEReaction
        self.kromeformat = self.reacformat.lower().strip()
        self.rate_string = None

        react_string = react_string.strip()
        if react_string != "" and react_string[0] != "#":
            kwords = self.kromeformat.split(",")

            for key, value in zip(kwords, react_string.split(",")):
                if value == "":
                    continue
                elif key == "idx":
                    self.idxfromfile = int(value)
                elif key == "r" and self._create_species(value):
                    self.reactants.append(self._create_species(value))
                elif key == "p" and self._create_species(value):
                    self.products.append(self._create_species(value))
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
