from textwrap import wrap, fill


def _stmwrap(text: str, width: int = 80, indent: int = 4):

    longindent = " " * indent
    shortindent = " " * (indent - 4)
    wrappedlist = wrap(text, width - indent, break_long_words=False)
    # wrappedstr = fill(text, width, subsequent_indent=longindent)
    wrappedstr = f"\n{longindent}".join(wrappedlist)
    wrappedstr = wrappedstr.replace(f"\n{longindent}}}", f"\n{shortindent}}}")

    return wrappedstr
