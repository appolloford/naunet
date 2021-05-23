from abc import ABC, abstractmethod
from collections import Counter
from enum import IntEnum
from ..species import Species

# If the reaction types have the same formalism, they share the same value in the enum class
class ReactionType(IntEnum):
    """
    The definition of reaction types. Types are categorized by phases
    (gas, grain, surface) and then divided to sub-groups
    """

    # 10x: common types of gas phase reactions
    GAS_TWOBODY = 100
    GAS_COSMICRAY = 101
    GAS_PHOTON = 102
    GAS_THREEBODY = 103

    # 11x: special types in KIDA
    GAS_KIDA_IP1 = 110
    GAS_KIDA_IP2 = 111

    # 12x: special types in UMIST
    GAS_UMIST_CRPHOT = 120

    # 13x: special types from Walsh+2015
    GAS_LEEDS_XRAY = 130
    GAS_LEEDS_RECOM = 131
    GAS_LEEDS_ECAPTURE = 132

    # 20x: common types of gas-grain interaction
    GRAIN_FREEZE = 200
    GRAIN_DESORPT_THERMAL = 201
    GRAIN_DESORPT_COSMICRAY = 202
    GRAIN_DESORPT_PHOTON = 203

    # 21x: special desorption from Walsh+2015
    GRAIN_DESORPT_REACTIVE = 210

    # 22x: special desoprtion from UCLCHEM
    GRAIN_DESORPT_H2 = 220

    # 30x: common types of surface reactions
    SURFACE_TWOBODY = 300
    SURFACE_COSMICRAY = 301
    SURFACE_PHOTON = 302

    UNKNOWN = 999


class Reaction(ABC):
    """
    Interface of reaction
    """

    consts = {}
    globs = {}
    vars = {}
    user_var = []

    def __init__(self, react_string, *args, **kwargs) -> None:
        self.reactants = []
        self.products = []
        self.temp_min = 1.0
        self.temp_max = -1.0
        self.reaction_type = ReactionType.UNKNOWN
        self.database = None

    def __str__(self) -> str:
        verbose = "{} -> {}".format(
            " + ".join(x.name for x in sorted(self.reactants)),
            " + ".join(x.name for x in sorted(self.products)),
        )
        return verbose

    def __eq__(self, o: object) -> bool:
        if not isinstance(o, Reaction):
            raise ValueError(f"{o} is not a reaction")

        # ignore reaction_type if there is no this information e.g. krome
        return (
            self.rpeq(o)
            and self.temp_min == o.temp_min
            and self.temp_max == o.temp_max
            and (
                self.reaction_type == o.reaction_type
                or self.reaction_type == ReactionType.UNKNOWN
                or o.reaction_type == ReactionType.UNKNOWN
            )
        )

    def __hash__(self) -> int:
        return hash(
            "_".join(str(x) for x in sorted(self.reactants) + sorted(self.products))
        )

    def __repr__(self) -> str:
        verbose = (
            (
                "{:16} -> {:32}, {:7.1f} < T < {:7.1f}, Type: {:25}, Database: {}".format(
                    " + ".join(x.name for x in self.reactants),
                    " + ".join(x.name for x in self.products),
                    self.temp_min,
                    self.temp_max,
                    self.reaction_type.name,
                    self.database,
                )
            )
            if len(self.reactants + self.products) > 0
            else " -> "
        )
        return verbose

    def _beautify(self, rate_string: str) -> str:
        """
        Beautify the reaction rate string

        :param rate_string: reaction rate string expression
        :type rate_string: str
        :return: beautified string
        :rtype: str
        """
        rate = (
            rate_string.replace("++", "+")
            .replace("--", "+")
            .replace("+-", "-")
            .replace("-+", "-")
        )
        return rate

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
        if species_name and species_name not in Species.known_pseudoelements():
            return Species(species_name)

    @classmethod
    def initialize(cls) -> None:
        """
        Interface. Change global settings / class attributes if needed
        """
        pass

    @classmethod
    def finalize(cls) -> None:
        """
        Interface. Reset global settings / class attributes if needed
        """
        pass

    @classmethod
    def preprocessing(cls, line: str) -> str:
        """
        Interface. Preprocess the input reaction string before initialize a reaction.
        Called in Network class

        :param line: reaction string
        :type line: str
        :return: proceeded string
        :rtype: str
        """
        return line

    def rpeq(self, o: object) -> bool:
        """
        Compare two reactions by their reactants and products.

        :param o: Another reaction
        :type o: object
        :return: True if two reactions have the same reactants and products
        :rtype: bool
        """
        return Counter(self.reactants) == Counter(o.reactants) and Counter(
            self.products
        ) == Counter(o.products)

    @abstractmethod
    def rate_func(self):
        """
        Abstract method. It should be implemented by child class.
        Return the rate function string.
        """
        raise NotImplementedError
