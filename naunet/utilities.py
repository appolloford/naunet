from __future__ import annotations
from collections import OrderedDict
from textwrap import wrap, fill
from .component import Component


def _collect_variable_items(
    complist: list[Component], var_type: str
) -> OrderedDict[str, str]:
    variables = OrderedDict()
    for comp in complist:
        var_dict = getattr(comp, var_type)
        for key, value in var_dict.items():
            variables[key] = value
    return variables.items()


def _prefix(text: str, pre: str) -> str:
    return "".join([pre, text])


def _suffix(text: str, suf: str) -> str:
    return "".join([text, suf])


def _stmwrap(text: str, width: int = 80, indent: int = 4):

    longindent = " " * indent
    shortindent = " " * (indent - 4)
    wrappedlist = wrap(text, width - indent, break_long_words=False)
    # wrappedstr = fill(text, width, subsequent_indent=longindent)
    wrappedstr = f"\n{longindent}".join(wrappedlist)
    wrappedstr = wrappedstr.replace(f"\n{longindent}}}", f"\n{shortindent}}}")

    return wrappedstr
