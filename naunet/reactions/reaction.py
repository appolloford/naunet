from __future__ import annotations
from collections import Counter
from ..dusts.dust import Dust
from ..species import Species
from .reactiontype import ReactionType


class Reaction:
    """Class of chemical reactions"""

    consts = {}
    varis = {
        "nH": None,
        "Tgas": None,
        "zeta": 1.3e-17,
        "Av": 1.0,
        "omega": 0.5,
    }
    locvars = []

    def __init__(
        self,
        reactants: list[Species] | list[str] = None,
        products: list[Species] | list[str] = None,
        temp_min: float = -1.0,
        temp_max: float = -1.0,
        alpha: float = 0.0,
        beta: float = 0.0,
        gamma: float = 0.0,
        reaction_type: ReactionType = ReactionType.UNKNOWN,
        format: str = "naunet",
        idxfromfile: int = -1,
        react_string: str = None,
    ) -> None:

        self.reactants = (
            [self._create_species(r) for r in reactants if self._create_species(r)]
            if reactants
            else []
        )
        self.products = (
            [self._create_species(p) for p in products if self._create_species(p)]
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

        self.react_string = react_string
        # if react_string:
        self._parse_string(react_string)

    def __contains__(self, spec: Species) -> bool:
        if not isinstance(spec, Species):
            return NotImplemented
        return spec in self.reactants or spec in self.products

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

        def fill(orig: list, nitem: int, dummy) -> list:
            return orig + [dummy] * (nitem - len(orig))

        if not format:
            verbose = str(self)

        elif format == "short":
            verbose = "{} -> {}".format(
                " + ".join(x.name for x in sorted(self.reactants)),
                " + ".join(x.name for x in sorted(self.products)),
            )

        elif format == "naunet":
            rnames = fill([f"{x:>12}" for x in sorted(self.reactants)], 3, f"{'':>12}")
            pnames = fill([f"{x:>12}" for x in sorted(self.products)], 5, f"{'':>12}")
            verbose = ",".join(
                [
                    f"{self.idxfromfile:<5}",
                    f",".join(rnames),
                    f",".join(pnames),
                    f"{self.alpha:10.3e}",
                    f"{self.beta:10.3e}",
                    f"{self.gamma:10.3e}",
                    f"{self.temp_min:9.2f}",
                    f"{self.temp_max:9.2f}",
                    f"{self.reaction_type:>4}",
                    f"{self.format:>8}",
                ]
            )

        elif format == "kida":
            rnames = fill([f"{x:<11}" for x in sorted(self.reactants)], 3, f"{'':>11}")
            pnames = fill([f"{x:<11}" for x in sorted(self.products)], 5, f"{'':>11}")
            verbose = " ".join(
                [
                    f"".join(rnames),
                    f"".join(pnames),
                    f"{self.alpha:+9.3e}" if self.alpha < 0 else f" {self.alpha:9.3e}",
                    f"{self.beta:+9.3e}" if self.beta < 0 else f" {self.beta:9.3e}",
                    f"{self.gamma:+9.3e}" if self.gamma < 0 else f" {self.gamma:9.3e}",
                    f"xxxxxxxx",
                    f"xxxxxxxx",
                    f"xxxx  x",
                    f"{int(self.temp_min):>6d}",
                    f"{int(self.temp_max):>6d}",
                    f"xx",
                    f"{self.idxfromfile:>5d}",
                    f"x  x",
                ]
            )

        elif format == "krome":
            rnames = fill([f"{x}" for x in sorted(self.reactants)], 3, "")
            pnames = fill([f"{x:<11}" for x in sorted(self.products)], 5, "")
            verbose = ",".join(
                [
                    f"{self.idxfromfile}",
                    f",".join(rnames),
                    f",".join(pnames),
                    f"{self.temp_min:.2f}",
                    f"{self.temp_max:.2f}",
                ]
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

    def _parse_string(self, react_string: str) -> None:

        if not react_string:
            return

        idx, *rps, a, b, c, lt, ut, rtype, format = react_string.split(",")
        self.reactants = [
            self._create_species(r.strip())
            for r in rps[0:3]
            if self._create_species(r.strip())
        ]
        self.products = [
            self._create_species(p.strip())
            for p in rps[3:8]
            if self._create_species(p.strip())
        ]

        self.alpha = float(a)
        self.beta = float(b)
        self.gamma = float(c)
        self.temp_min = float(lt)
        self.temp_max = float(ut)
        self.idxfromfile = int(idx)
        self.reaction_type = ReactionType(int(rtype))
        self.format = format

    def _create_species(self, species_name: Species | str) -> Species:
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

    @property
    def is_empty(self) -> bool:
        """
        Check whether the reaction if empty

        Returns:
            bool: True if there is neither reactants nor products
        """
        return len(self.reactants) == 0 and len(self.products) == 0

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
            rate = " * ".join(
                s
                for s in [
                    f"{a}",
                    f"pow(Tgas/300.0, {b})" if b else "",
                    f"exp(-{c}/Tgas)" if c else "",
                ]
                if s
            )

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

        elif rtype == ReactionType.DUMMY:
            rate = "0.0"

        else:
            raise RuntimeError(f"Unknown reaction type {self.reaction_type}")

        return rate
