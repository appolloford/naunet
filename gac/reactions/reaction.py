from abc import ABC, abstractmethod
from enum import Enum
from ..species import Species
from ..settings import pseudo_element_list

# If the reaction types have the same formalism, they share the same value in the enum class
class ReactionType(Enum):
    COSMIC_RAY = 1
    PHOTON = 2
    TWO_BODY = 3

    KIDA_CR = 1  # Cosmic-ray ionization
    KIDA_PD = 2  # Photo-dissociation (Draine)
    KIDA_MA = 3  # Modified Arrhenius
    KIDA_IP1 = 4  # ionpol1
    KIDA_IP2 = 5  # ionpol2
    KIDA_TB = 6  # Three-body

    UMIST_AD = 3  # Associative Detachment
    UMIST_CD = 3  # Collisional Dissociation
    UMIST_CE = 3  # Charge Exchange
    UMIST_CP = 1  # Cosmic-Ray Proton (CRP)
    UMIST_CR = 7  # Cosmic-Ray Photon (CRPHOT)
    UMIST_DR = 3  # Dissociative Recombination
    UMIST_IN = 3  # Ion-Nuetral
    UMIST_MN = 3  # Mutual Neutralisation
    UMIST_NN = 3  # Nuetral-Neutral
    UMIST_PH = 2  # Photoprocess
    UMIST_RA = 3  # Radiative Association
    UMIST_REA = 3  # Radiative Electron Attachment
    UMIST_RR = 3  # Radiative Recombination


class Reaction(ABC):
    def __init__(self, react_string) -> None:
        self.reactants = []
        self.products = []
        self.temp_min = 1.0
        self.temp_max = -1.0
        self.reaction_type = None
        self.database = None

        self._parse_string(react_string)

    def __str__(self) -> str:
        verbose = "{} -> {}".format(
            " + ".join(self.reactants), " + ".join(self.products)
        )
        return super().__str__() + verbose

    def __eq__(self, o: object) -> bool:
        return (
            self.rpeq(o)
            and self.reaction_type == o.reaction_type
            and self.temp_min == o.temp_min
            and self.temp_max == o.temp_max
        )

    @abstractmethod
    def _parse_string(self, react_string) -> None:
        """
        Abstract method. It should be implemented by child class.
        Called by __init__. Parsing the input string to set the
        reactants, product, reaction rate function, etc.
        """
        raise NotImplementedError

    def create_species(self, species_name: str) -> object:
        """
        Create a Species instance if the name is not a pseudo element (e.g. CR, CRPHOT), else return None

        :param species_name: the name of the species
        :type species_name: str
        :return: Species object
        :rtype: object
        """
        if species_name not in pseudo_element_list:
            return Species(species_name)

    def rpeq(self, o: object) -> bool:
        """
        Compare two reactions by their reactants and products.

        :param o: Another reaction
        :type o: object
        :return: True if two reactions have the same reactants and products
        :rtype: bool
        """
        req = sorted(self.reactants, key=lambda p: p.name) == sorted(
            o.reactants, key=lambda p: p.name
        )
        peq = sorted(self.products, key=lambda p: p.name) == sorted(
            o.products, key=lambda p: p.name
        )
        return req and peq

    @abstractmethod
    def rate_func(self) -> str:
        """
        Abstract method. It should be implemented by child class.
        Return the rate function string.
        """
        raise NotImplementedError
