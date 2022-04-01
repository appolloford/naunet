from __future__ import annotations
from collections import Counter
from enum import IntEnum
from ..dusts.dust import Dust
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

    # 31x: surface diffusion
    SURFACE_DIFFUSION = 310

    UNKNOWN = 999


class Reaction:
    """Class of chemical reactions"""

    consts = {}
    varis = {
        "nH": None,
        "Tgas": None,
        "zeta": 1.3e-17,
    }
    locvars = []

    def __init__(
        self,
        reactants: list[Species] | list[str] = None,
        products: list[Species] | list[str] = None,
        temp_min: float = 1.0,
        temp_max: float = -1.0,
        alpha: float = 0.0,
        beta: float = 0.0,
        gamma: float = 0.0,
        reaction_type: ReactionType = ReactionType.UNKNOWN,
        format: str = None,
        idxfromfile: int = -1,
    ) -> None:

        self.reactants = (
            [self.create_species(r) for r in reactants if self.create_species(r)]
            if reactants
            else []
        )
        self.products = (
            [self.create_species(p) for p in products if self.create_species(p)]
            if products
            else []
        )
        self.temp_min = temp_min
        self.temp_max = temp_max
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.reaction_type = reaction_type
        self.format = format
        self.idxfromfile = idxfromfile

    def __str__(self) -> str:
        verbose = (
            (
                "{:16} -> {:32}, {:7.1f} < T < {:7.1f}, Type: {:25}, Format: {}, Index: {}".format(
                    " + ".join(x.name for x in self.reactants),
                    " + ".join(x.name for x in self.products),
                    self.temp_min,
                    self.temp_max,
                    self.reaction_type.name,
                    self.format,
                    self.idxfromfile,
                )
            )
            if len(self.reactants + self.products) > 0
            else " -> "
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

    def __format__(self, format: str) -> str:

        verbose = None

        if not format:
            verbose = str(self)

        elif format == "short":
            verbose = "{} -> {}".format(
                " + ".join(x.name for x in sorted(self.reactants)),
                " + ".join(x.name for x in sorted(self.products)),
            )

        return verbose

    def __hash__(self) -> int:
        return hash(
            "_".join(str(x) for x in sorted(self.reactants) + sorted(self.products))
        )

    def __repr__(self) -> str:
        params = ", ".join(
            [
                f"{self.reactants}",
                f"{self.products}",
                f"{self.temp_min}",
                f"{self.temp_max}",
                f"{self.alpha}",
                f"{self.beta}",
                f"{self.gamma}",
                f"ReactionType.{self.reaction_type.name}",
                f"'{self.format}'",
                f"{self.idxfromfile}",
            ]
        )

        verbose = f"Reaction({params})"
        return verbose

    def _beautify(self, rate_string: str) -> str:
        """
        Beautify the reaction rate string

        Args:
            rate_string (str): reaction rate expression

        Returns:
            str: beautified reaction rate expression
        """

        rate = (
            rate_string.replace("++", "+")
            .replace("--", "+")
            .replace("+-", "-")
            .replace("-+", "-")
        )
        return rate

    def create_species(self, species_name: Species | str) -> Species:
        """
        Create a Species instance if the name is not a pseudo element
        (e.g. CR, CRPHOT), else return None

        Args:
            species_name (Species | str): name of species

        Returns:
            Species: Species instance of the input name
        """

        if isinstance(species_name, Species):
            return species_name

        if species_name and species_name not in Species.known_pseudoelements():
            return Species(species_name)

    @classmethod
    def initialize(cls) -> None:
        """
        Change settings / class attributes if needed
        """

        pass

    @classmethod
    def finalize(cls) -> None:
        """
        Reset settings / class attributes if needed
        """

        pass

    @classmethod
    def preprocessing(cls, line: str) -> str:
        """
        Preprocess the input reaction string before initialize a reaction.
        Called in Network class to deal with input with special meanings.

        Args:
            line (str): input string of a reaction

        Returns:
            str: proceeded input string
        """

        return line

    def rpeq(self, o: object) -> bool:
        """
        Compare two reactions by their reactants and products.

        Args:
            o (object): other reaction instance

        Returns:
            bool: True if the reactants and products are the same in two
            reactions. Otherwise False
        """

        return Counter(self.reactants) == Counter(o.reactants) and Counter(
            self.products
        ) == Counter(o.products)

    def rateexpr(self, dust: Dust = None) -> str:
        """Returns the reaction rate expression in C language

        Args:
            dust (Dust, optional): Dust model to be used. Defaults to None.

        Raises:
            RuntimeError: if the reaction type is unknown

        Returns:
            str: the reaction rate expression in C language
        """
        a = self.alpha
        b = self.beta
        c = self.gamma

        rtype = self.reaction_type

        # two-body gas-phase reaction
        if rtype == ReactionType.GAS_TWOBODY:
            rate = f"{a} * pow(Tgas/300.0, {b}) * exp(-{c}/Tgas)"

        elif rtype == ReactionType.GAS_COSMICRAY:
            rate = f"{a} * zeta"

        elif rtype == ReactionType.GAS_PHOTON:
            rate = f"{a} * exp(-{c}*Av)"

        elif rtype == ReactionType.GAS_KIDA_IP1:
            rate = f"{a} * {b} * (0.62 + 0.4767*{c}*sqrt(300.0/Tgas))"

        elif rtype == ReactionType.GAS_KIDA_IP2:
            rate = f"{a} * {b} * (1 + 0.0967*{c}*sqrt(300.0/Tgas) + {c}*{c}*(300.0/Tgas)/10.526)"

        elif rtype == ReactionType.GAS_UMIST_CRPHOT:
            rate = f"{a} * pow(Tgas/300.0, {b}) * {c} / (1-omega)"

        else:
            raise RuntimeError(f"Unknown reaction type {self.reaction_type}")

        return rate
