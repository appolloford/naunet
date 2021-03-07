import sys


def conti_confirm(question: str, default: bool = False) -> None:
    answer = input(f"{question} [Y/n]" if default else f"{question} [y/N]")

    if answer.lower().startswith("n"):
        sys.exit()

    elif not answer and default == False:
        sys.exit()
