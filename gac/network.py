import logging
import sympy as sym
from operator import mul
from functools import reduce
from .species import Species
from .reactions.reaction import Reaction
from .reactions.kidareaction import KIDAReaction

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()


def reaction_factory(react_string: str, database: str) -> Reaction:
    if database == "kida":
        return KIDAReaction(react_string)


class Network:
    def __init__(self, fname: str = None, database: str = None) -> None:
        self.reaction_list = []
        self.reactants_in_network = []
        self.products_in_network = []

        if fname and database:
            self.add_reaction_from_file(fname, database)

    def _add_reaction(self, react_string: str, database: str) -> list:
        reaction = reaction_factory(react_string, database)
        self.reaction_list.append(reaction)
        new_reactants = [
            spec for spec in reaction.reactants if spec not in self.reactants_in_network
        ]
        new_products = [
            spec for spec in reaction.products if spec not in self.products_in_network
        ]
        self.reactants_in_network.extend(new_reactants)
        self.products_in_network.extend(new_products)
        if len(self.reaction_list) % 100 == 0:
            print("Processing: {} reactions...".format(len(self.reaction_list)))
        return new_reactants + new_products

    def add_reaction(self, react_string: str, database: str) -> None:
        new_species = self._add_reaction(react_string, database)
        logger.info("New species are added: {}".format(new_species))

    def add_reaction_from_file(self, filename: str, database: str) -> None:
        if not database:
            logger.warning(
                'Try to read in file but database is not assigned. Try again by "add_reaction_from_file"'
            )
        with open(filename, "r") as networkfile:
            new_species = []
            for line in networkfile.readlines():
                new_species.extend(self._add_reaction(line, database))
            for x in new_species:
                print(x)

    def check_duplicate_reaction(self):
        dupes = [
            x for n, x in enumerate(self.reaction_list) if x in self.reaction_list[:n]
        ]
        logger.info(
            "The following reactions are duplicate: {}".format("\n".join(dupes))
        )
        print("The following reactions are duplicate:", dupes)
        return dupes

    def check_source_sink(self):
        source = [
            x for x in self.reactants_in_network if x not in self.products_in_network
        ]
        sink = [
            x for x in self.products_in_network if x not in self.reactants_in_network
        ]
        if source == [] and sink == []:
            print("Found no source or sink")
        elif source != []:
            print("Found sources: ", source)
        elif sink != []:
            print("Found sinks: ", sink)

    def to_ccode(self):
        reacprod = self.reactants_in_network + self.products_in_network
        species_in_network = [
            x for n, x in enumerate(reacprod) if x not in reacprod[:n]
        ]
        nspec = len(species_in_network)
        ydot = sym.MatrixSymbol("ydot", nspec, 1)
        y = sym.MatrixSymbol("y", nspec, 1)
        rhs = [0.0] * nspec
        for react in self.reaction_list:
            rate = react.rate_func()
            ridx = [species_in_network.index(r) for r in react.reactants]
            pidx = [species_in_network.index(p) for p in react.products]
            rsym = [y[idx] for idx in ridx]
            # psym = [y[idx] for idx in pidx]
            rsym_product = reduce(mul, rsym)
            # psym_product = reduce(mul, psym)
            for idx in ridx:
                rhs[idx] -= rate * rsym_product
            for idx in pidx:
                rhs[idx] += rate * rsym_product

        for left, right in zip(ydot, rhs):
            print(sym.ccode(right, assign_to=left))
