import logging
from sympy import symbols, MatrixSymbol, Idx, ccode
from sympy.codegen.ast import Assignment
from sympy.utilities.codegen import codegen
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
        self.net_species = []
        self.nspecies = 0
        self.info_updated = False

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
        self.info_updated = False

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

        self.info_updated = False

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

    def collect_infos(self):
        if self.info_updated:
            return

        reacprod = self.reactants_in_network + self.products_in_network
        self.net_species = [x for n, x in enumerate(reacprod) if x not in reacprod[:n]]
        self.nspecies = len(self.net_species)
        logging.info(
            "{} species in the network: {}".format(
                self.nspecies, ", ".join([x.name for x in self.net_species])
            )
        )
        self.info_updated = True

    def fex(self):
        rhs = self.rhs()
        ydot = MatrixSymbol("ydot", self.nspecies, 1)

        fex = [Assignment(l, r) for l, r in zip(ydot, rhs)]
        return fex

    def fex_to_ccode(
        self, to_file=False, prefix="fex", ext="c", header=False, h_ext="h"
    ):

        tab = "    "

        if to_file:

            cfunc = (
                "static int fex(realtype t, N_Vector y, N_Vector ydot, void *user_data)"
            )

            with open(f"{prefix}.{ext}", "w") as cfile:
                if not header:
                    cfile.write("#include <cvode/cvode.h>\n")
                    cfile.write("#include <nvector/nvector_serial.h>\n\n")
                else:
                    cfile.write(f'#include "{prefix}.{h_ext}"\n\n')
                cfile.write(f"{cfunc}")
                cfile.write("{\n\n")
                cfile.write(f"{tab}UserData *u_data = (UserData*) user_data;\n")
                cfile.write(f"{tab}realtype Tgas = u_data->Tgas;\n\n")
                cfile.write("\n".join([tab + ccode(f) for f in self.fex()]))
                cfile.write("\n\n}")

            if header:
                with open(f"{prefix}.{h_ext}", "w") as hfile:
                    hfile.write("#ifndef __FEX_H__\n")
                    hfile.write("#define __FEX_H__\n\n")
                    hfile.write("#include <cvode/cvode.h>\n")
                    hfile.write("#include <nvector/nvector_serial.h>\n\n")
                    hfile.write(f"{cfunc};\n\n")
                    hfile.write("#endif\n")
        else:
            print("\n".join([ccode(f) for f in self.fex()]))

    def rhs(self):
        if not self.info_updated:
            self.collect_infos()

        y = MatrixSymbol("y", self.nspecies, 1)
        rhs = [0.0] * self.nspecies

        for react in self.reaction_list:
            rate = react.rate_func()
            ridx = [self.net_species.index(r) for r in react.reactants]
            pidx = [self.net_species.index(p) for p in react.products]
            rsym = [y[idx] for idx in ridx]
            # psym = [y[idx] for idx in pidx]
            rsym_product = reduce(mul, rsym)
            # psym_product = reduce(mul, psym)
            for idx in ridx:
                rhs[idx] -= rate * rsym_product
            for idx in pidx:
                rhs[idx] += rate * rsym_product

        return rhs
