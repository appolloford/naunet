from lark import Lark, Transformer


class ExpressionConverter:
    class Expression(Transformer):
        expression = lambda self, e: " ".join(e)
        multiply = lambda self, m: "".join(m)
        power = lambda self, p: "".join(p)
        func = lambda self, f: "".join(f)
        variable = lambda self, v: "".join(v)
        listvar = lambda self, l: "".join(l)
        index = lambda self, i: f"IDX{''.join(i)}"
        atom = lambda self, a: "".join(a)
        COMMA = lambda self, c: ", "
        TIMES = lambda self, t: " * "

        def scientific(self, s):
            (s,) = s
            return s.value

    class FExpression(Expression):
        power = (
            lambda self, p: "".join(p)
            .replace("pow(", "(")
            .replace(")", "")
            .replace(", ", ")**")
        )
        listvar = (
            lambda self, l: "".join(l)
            .replace("[", "(")
            .replace("]", ")")
            .replace("y", "n")
        )
        index = lambda self, i: f"idx{''.join(i)}"

    class CExpression(Expression):
        power = lambda self, p: f"pow({''.join(p).replace('**', ', ')})"
        listvar = (
            lambda self, l: "".join(l)
            .replace("(", "[")
            .replace(")", "]")
            .replace("n", "y")
        )

        # def index(self, i):
        #     (i,) = i
        #     return i.value

    fgrammar = r"""
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

    cgrammar = r"""
        expression: multiply ((PLUS | MINUS) multiply)*
        multiply: atom ((TIMES | DIV) atom)*
        power: POW LPAREN atom COMMA atom RPAREN
        listvar: variable LSQBKT index RSQBKT
        func: variable LPAREN expression (COMMA expression)* RPAREN
        variable: WORD (WORD | NUMBER | UNDER)*
        index: "IDX" UNDER WORD (WORD| NUMBER | UNDER)*
        scientific: NUMBER ((E1 | E2) SIGN? NUMBER)?
        atom: scientific
            | power
            | variable
            | listvar
            | func
            | LPAREN expression RPAREN
            | expression
        PLUS: "+"
        MINUS: "-"
        TIMES: "*"
        DIV: "/"
        POW: "pow"
        LPAREN: "("
        RPAREN: ")"
        LSQBKT: "["
        RSQBKT: "]"
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

    grammar = {
        "fortran": fgrammar,
        "c": cgrammar,
    }

    transformer = {
        "fortran": FExpression,
        "c": CExpression,
    }

    def __init__(self, lang: str) -> None:
        self._lang = lang.lower()
        grammar = self.grammar.get(self._lang)
        self._parser = Lark(grammar, start="expression")
        self._expr = None

    def read(self, expression: str) -> None:
        self._expr = self._parser.parse(expression)

    def __str__(self) -> str:
        transformer = self.transformer.get("c")
        expr = transformer().transform(self._expr)
        return expr

    def __format__(self, format: str) -> str:
        format = format.lower() if format else self._lang
        transformer = self.Expression
        if format != self._lang:
            transformer = self.transformer.get(format)
        expr = transformer().transform(self._expr)
        return expr
