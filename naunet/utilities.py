from textwrap import wrap, fill


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
