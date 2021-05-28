from abc import ABC, abstractmethod


class Dust(ABC):

    consts = {}
    globs = {}
    varis = {}
    user_var = []

    def __init__(self, *args, **kwargs) -> None:
        super().__init__()
