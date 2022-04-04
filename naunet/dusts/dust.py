from abc import ABC, abstractmethod


class Dust(ABC):

    consts = {}
    varis = {}
    locvars = []

    def __init__(self, *args, **kwargs) -> None:
        self.model = "none"
        super().__init__()

    @abstractmethod
    def rate_depletion(self, *args, **kwargs) -> str:
        raise NotImplementedError

    @abstractmethod
    def rate_desorption(self, *args, **kwargs) -> str:
        raise NotImplementedError

    @abstractmethod
    def rate_electroncapture(self, *args, **kwargs) -> str:
        raise NotImplementedError

    @abstractmethod
    def rate_recombination(self, *args, **kwargs) -> str:
        raise NotImplementedError

    @abstractmethod
    def rate_surface1(self, *args, **kwargs) -> str:
        """
        One-species surface reaction rate
        """
        raise NotImplementedError

    @abstractmethod
    def rate_surface2(self, *args, **kwargs) -> str:
        """
        Two-species surface reaction rate
        """
        raise NotImplementedError
