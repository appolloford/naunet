import logging
import re
from lark import Lark, Transformer
from .. import settings
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
            .replace("n", settings.ode_symbols["ode_vector"])
        )
        index = lambda self, i: "idx" + "".join(i)
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


class Expression(Transformer):
    def powexpression(self, p):
        return f"pow({''.join(p).replace('**', ', ')})"

    def expression(self, e):
        return " ".join(e)

    def scientific(self, s):
        (s,) = s
        return s.value

    def atom(self, a):
        return "".join(a)

    def func(self, f):
        return "".join(f)

    def WORD(self, w):
        return w.value

    def generalexpression(self, g):
        return "".join(g)

    def variable(self, f):
        return "".join(f)


class KROMEReaction(Reaction):

    reacformat = "idx,r,r,r,p,p,p,p,tmin,tmax,rate"
    common = []
    var = []

    def __init__(self, react_string, *args, **kwargs) -> None:
        super().__init__(react_string)

        self.database = "krome"
        self.alpha = 0.0
        self.beta = 0.0
        self.gamma = 0.0
        self.kromeformat = self.reacformat.lower().strip()
        self.rate_string = None

        self._parse_string(react_string)

    @classmethod
    def preprocessing(cls, line: str) -> str:
        if line.startswith(("#", "//")):
            return ""
        elif line.startswith("@format:"):
            KROMEReaction.reacformat = line.replace("@format:", "")
            return ""
        elif line.startswith("@var"):
            KROMEReaction.var.extend(line.replace("@var:", "").split(","))
        elif line.startswith("@common:"):
            KROMEReaction.common.extend(line.replace("@common:", "").split(","))
            return ""
        else:
            return line.strip()

    def rate_func(self):

        # grammar = r"""
        #     generalexpression: expression ((PLUS | MINUS) expression)*
        #     expression: atom ((TIMES | DIV) atom)*
        #     atom: powexpression
        #         | scientific
        #         | LPAREN expression RPAREN
        #         | LPAREN generalexpression RPAREN
        #         | WORD
        #         | func
        #         | variable
        #     powexpression: (atom POW atom)
        #     func: variable LPAREN generalexpression (COMMA generalexpression)* RPAREN
        #     PLUS: "+"
        #     MINUS: "-"
        #     TIMES: "*"
        #     DIV: "/"
        #     POW: "**"
        #     LPAREN: "("
        #     RPAREN: ")"
        #     COMMA: ","
        #     variable: WORD (WORD | NUMBER)*
        #     scientific: SCIENTIFIC_NUMBER
        #     SCIENTIFIC_NUMBER: NUMBER ((E1 | E2) SIGN? NUMBER)?
        #     E1: "E"
        #     E2: "e"
        #     SIGN: "+" | "-"
        #     %import common.WORD             -> WORD
        #     %import common.SIGNED_NUMBER    -> NUMBER
        #     %import common.WS
        #     %ignore WS
        # """

        # rate_parser = Lark(grammar, start="generalexpression")

        # print(type(self.rate_string), self.rate_string)
        rate = re.sub(r"(\d\.?)d(\-?\d)", r"\1e\2", self.rate_string)
        # rate = rate_parser.parse(rate)
        # rate = Expression().transform(rate)
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
