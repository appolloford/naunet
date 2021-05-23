import logging
import re
from lark import Lark, Transformer
from ..species import Species
from .reaction import Reaction, ReactionType


class FtoCCoverter:
    class CExpression(Transformer):
        expression = lambda self, e: " ".join(e)
        multiply = lambda self, m: "".join(m)
        power = lambda self, p: f"pow({''.join(p).replace('**', ', ')})"
        func = lambda self, f: "".join(f)
        variable = lambda self, v: "".join(v)
        listvar = (
            lambda self, l: "".join(l)
            .replace("(", "[")
            .replace(")", "]")
            .replace("n", "y")
        )
        index = lambda self, i: "IDX" + "".join(i)
        atom = lambda self, a: "".join(a)

        def scientific(self, s):
            (s,) = s
            return s.value

        # def index(self, i):
        #     (i,) = i
        #     return i.value

    transformer = CExpression()

    grammar = r"""
        expression: multiply ((PLUS | MINUS) multiply)*
        multiply: atom ((TIMES | DIV) atom)*
        power: (atom POW atom)
        func: variable LPAREN expression (COMMA expression)* RPAREN
        variable: WORD (WORD | NUMBER | UNDER)*
        listvar: variable LPAREN index RPAREN
        index: "idx" UNDER WORD (WORD| NUMBER | UNDER)*
        scientific: NUMBER ((E1 | E2) SIGN? NUMBER)?
        atom: scientific
            | power
            | variable
            | listvar
            | func
            | LPAREN expression RPAREN
        PLUS: "+"
        MINUS: "-"
        TIMES: "*"
        DIV: "/"
        POW: "**"
        LPAREN: "("
        RPAREN: ")"
        COMMA: ","
        E1: "E"
        E2: "e"
        SIGN: "+" | "-"
        UNDER: "_"
        %import common.WORD             -> WORD
        %import common.SIGNED_NUMBER    -> NUMBER
        %import common.WS
        %ignore WS
    """

    parser = Lark(grammar, start="expression")

    def __init__(self) -> None:
        pass

    @classmethod
    def convert(cls, expression: str) -> str:
        cexpr = cls.parser.parse(expression)
        cexpr = cls.transformer.transform(cexpr)
        return cexpr


class KROMEReaction(Reaction):

    reacformat = "idx,r,r,r,p,p,p,p,tmin,tmax,rate"
    vars = {
        "Hnuclei": "nH",
        "Temperature": "Tgas",
    }
    user_var = []

    def __init__(self, react_string, *args, **kwargs) -> None:
        super().__init__(react_string)

        self.database = "KROME"
        self.alpha = 0.0
        self.beta = 0.0
        self.gamma = 0.0
        self.kromeformat = self.reacformat.lower().strip()
        self.rate_string = None

        self._parse_string(react_string)

    @classmethod
    def finalize(cls) -> None:
        # restore the default settings after completing a file
        print("krome finalize is called")
        cls.reacformat = "idx,r,r,r,p,p,p,p,tmin,tmax,rate"
        cls.vars = {
            "Hnuclei": "nH",
            "Temperature": "Tgas",
        }
        cls.user_var = []

    @classmethod
    def preprocessing(cls, line: str) -> str:
        if line.startswith(("#", "//")):
            return ""
        elif line.startswith("@format:"):
            cls.reacformat = line.replace("@format:", "")
            return ""
        elif line.startswith("@var"):
            if "Hnuclei" not in line:
                cls.user_var.append(line.replace("@var:", "").strip())
            return ""
        elif line.startswith("@common:"):
            commonlist = line.replace("@common:", "").strip().split(",")
            cls.vars.update(zip(commonlist, commonlist))
            return ""
        else:
            return line.strip()

    def rate_func(self):

        rate = re.sub(r"(\d\.?)d(\-?\d)", r"\1e\2", self.rate_string)
        rate = re.sub(r"(idx_.?)p", r"\1II", rate)
        rate = re.sub(r"(idx_.?)m", r"\1M", rate)
        rate = re.sub(r"(idx_.?)\)", r"\1I)", rate)
        rate = rate.replace("Hnuclei", "nH")
        rate = FtoCCoverter.convert(rate)
        return rate

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
