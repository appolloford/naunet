import logging
from sympy import symbols, Symbol, Function, MatrixSymbol, Idx, ccode
from sympy.codegen.ast import Assignment, CodeBlock, Declaration, Variable, real
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
        self.reactants_in_network = set()
        self.products_in_network = set()
        self.net_species = []
        self.nspecies = 0
        self.info_updated = False

        if fname and database:
            self.add_reaction_from_file(fname, database)

    def _add_reaction(self, react_string: str, database: str) -> list:
        reaction = reaction_factory(react_string, database)
        self.reaction_list.append(reaction)
        new_reactants = reaction.reactants.difference(self.reactants_in_network)
        new_products = reaction.products.difference(self.products_in_network)
        self.reactants_in_network.update(new_reactants)
        self.products_in_network.update(new_products)
        if len(self.reaction_list) % 100 == 0:
            print("Processing: {} reactions...".format(len(self.reaction_list)))
        return new_reactants | new_products

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
            new_species = set()
            for line in networkfile.readlines():
                new_species.update(self._add_reaction(line, database))
            for x in new_species:
                print(x)

        self.info_updated = False

    def check_duplicate_reaction(self, full_check: bool = True):

        seen = {}
        dupes = []

        check_list = (
            self.reaction_list
            if full_check
            else [str(react) for react in self.reaction_list]
        )

        for x in check_list:
            if x not in seen:
                seen[x] = 1
            else:
                if seen[x] == 1:
                    dupes.append(x)
                seen[x] += 1

        logger.info(
            "The following reactions are duplicate:\n{}".format(
                "\n".join([repr(x) for x in dupes])
            )
        )
        print(
            "The following reactions are duplicate:\n{}".format(
                "\n".join([repr(x) for x in dupes])
            )
        )
        return dupes

    def check_source_sink(self):
        source = self.reactants_in_network.difference(self.products_in_network)
        sink = self.products_in_network.difference(self.reactants_in_network)
        if len(source) == 0 and len(sink) == 0:
            print("Found no source or sink")
        elif len(source) != 0:
            print("Found sources: ", source)
        elif len(sink) != 0:
            print("Found sinks: ", sink)

    def collect_infos(self):
        if self.info_updated:
            return

        self.net_species = list(self.reactants_in_network | self.products_in_network)
        self.nspecies = len(self.net_species)
        logging.info(
            "{} species in the network: {}".format(
                self.nspecies, ", ".join([x.name for x in self.net_species])
            )
        )
        self.info_updated = True

    # TODO: complete the function
    def constants_header(self):
        if not self.info_updated:
            self.collect_infos()

        with open("test/test_output/constants.h", "w") as hfile:
            hfile.write("#ifndef __CONSTANT_H__\n")
            hfile.write("#define __CONSTANT_H__\n\n")
            hfile.write(f"#define NSPECIES {self.nspecies}\n\n")
            hfile.write(
                "\n".join(
                    f"#define IDX_{x.alias} {i}" for i, x in enumerate(self.net_species)
                )
            )
            hfile.write("\n\n")
            hfile.write("#endif\n")

    def fex(self, symbol_form=False):
        rate_declare, rate_assign, rhs = self.rhs(symbol_form)

        # TODO: Return a pure sympy object list when it is possible, make it can be handle by ccode()
        # * Sympy doesn't generate semicolons after declaration, neither support the if statement.
        # * Therefore, fex() return the declaration and assignment with the final function.
        # * fex_to_code() will do further processing to generate the correct code.
        # fex_result = rate_declare
        # fex_result += rate_assign

        ydot = MatrixSymbol("ydot", self.nspecies, 1)
        fex_result = [Assignment(l, r) for l, r in zip(ydot, rhs)]
        return rate_declare, rate_assign, fex_result

    def fex_to_ccode(
        self,
        *args,
        to_file=False,
        prefix="",
        fname="fex",
        ext="c",
        header=False,
        h_ext="h",
        **kwargs,
    ):

        tab = "    "

        if to_file:

            cfunc = "int fex(realtype t, N_Vector u, N_Vector u_dot, void *user_data)"
            rate_declare, rate_assign, fex = self.fex(*args, **kwargs)

            with open(f"{prefix + fname}.{ext}", "w") as cfile:
                if header:
                    cfile.write(f'#include "{fname}.{h_ext}"\n\n')
                    with open(f"{prefix + fname}.{h_ext}", "w") as hfile:
                        hfile.write("#ifndef __JAC_H__\n")
                        hfile.write("#define __JAC_H__\n\n")
                        hfile.write("#include <cvode/cvode.h>\n")
                        hfile.write("#include <nvector/nvector_serial.h>\n\n")
                        hfile.write(f"{cfunc};\n\n")
                        hfile.write("#endif\n")
                else:
                    cfile.write("#include <cvode/cvode.h>\n")
                    cfile.write("#include <nvector/nvector_serial.h>\n\n")

                cfile.write(
                    f"// Simple function that calculates the differential equation.\n"
                )
                cfile.write(f"{cfunc}")
                cfile.write("{\n\n")
                cfile.write(f"{tab}realtype *y = N_VGetArrayPointer(u);\n")
                cfile.write(f"{tab}realtype *ydot = N_VGetArrayPointer(u_dot);\n")
                cfile.write(f"{tab}UserData *u_data = (UserData*) user_data;\n")
                cfile.write(f"{tab}realtype Tgas = u_data->Tgas;\n\n")
                cfile.write(
                    "\n".join([f"{tab}{ccode(dec)} = 0.0;" for dec in rate_declare])
                )
                cfile.write("\n\n")
                cfile.write(
                    "\n".join(
                        [
                            f"{tab}if (Tgas>{react.temp_min} && Tgas<{react.temp_max}) {ccode(assign)}"
                            for react, assign in zip(self.reaction_list, rate_assign)
                        ]
                    )
                )
                cfile.write("\n\n")
                cfile.write("\n".join([tab + ccode(f) for f in fex]))
                # cfile.write(ccode(CodeBlock(*fex)))
                cfile.write("\n\n")
                cfile.write(f"{tab}return 0;")
                cfile.write("\n\n}")

        else:
            print("\n".join([ccode(f) for f in self.fex()]))

    def jac(self, symbol_form=False):

        y = MatrixSymbol("y", self.nspecies, 1)
        jac_result = [0.0] * self.nspecies * self.nspecies
        rate_declare = []
        rate_assign = []

        for r, react in enumerate(self.reaction_list):
            if symbol_form:
                rate_symbol = Symbol(f"rate_react{r}")
                rate_declare.append(Declaration(Variable(rate_symbol, type=real)))
                rate_assign.append(Assignment(rate_symbol, react.rate_func()))
            else:
                rate_symbol = react.rate_func()
            ridx = [self.net_species.index(r) for r in react.reactants]
            pidx = [self.net_species.index(p) for p in react.products]

            for idx in ridx:
                for ri in ridx:
                    rsym = [y[i] for i in ridx if i != ri]
                    rsym_product = reduce(mul, rsym)
                    jac_result[idx * self.nspecies + ri] -= rate_symbol * rsym_product
            for idx in pidx:
                for ri in ridx:
                    rsym = [y[i] for i in ridx if i != ri]
                    rsym_product = reduce(mul, rsym)
                    jac_result[idx * self.nspecies + ri] += rate_symbol * rsym_product

        return rate_declare, rate_assign, jac_result

    def jac_to_ccode(
        self,
        *args,
        to_file=False,
        prefix="",
        fname="jac",
        ext="c",
        header=False,
        h_ext="h",
        **kwargs,
    ):
        tab = "    "

        if to_file:

            cfunc = "int jac(realtype t, N_Vector u, N_Vector fu, SUNMatrix Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)"

            rate_declare, rate_assign, jac = self.jac(*args, **kwargs)

            with open(f"{prefix + fname}.{ext}", "w") as cfile:
                if header:
                    cfile.write(f'#include "{fname}.{h_ext}"\n\n')
                    with open(f"{prefix + fname}.{h_ext}", "w") as hfile:
                        hfile.write("#ifndef __FEX_H__\n")
                        hfile.write("#define __FEX_H__\n\n")
                        hfile.write("#include <cvode/cvode.h>\n")
                        hfile.write("#include <nvector/nvector_serial.h>\n\n")
                        hfile.write(f"{cfunc};\n\n")
                        hfile.write("#endif\n")
                else:
                    cfile.write("#include <cvode/cvode.h>\n")
                    cfile.write("#include <nvector/nvector_serial.h>\n\n")

                cfile.write(f"{cfunc}")
                cfile.write("{\n\n")
                cfile.write(f"{tab}realtype *y = N_VGetArrayPointer(u);\n")
                cfile.write(f"{tab}UserData *u_data = (UserData*) user_data;\n")
                cfile.write(f"{tab}realtype Tgas = u_data->Tgas;\n\n")
                cfile.write(
                    "\n".join([f"{tab}{ccode(dec)} = 0.0;" for dec in rate_declare])
                )
                cfile.write("\n\n")
                cfile.write(
                    "\n".join(
                        [
                            f"{tab}if (Tgas>{react.temp_min} && Tgas<{react.temp_max}) {ccode(assign)}"
                            for react, assign in zip(self.reaction_list, rate_assign)
                        ]
                    )
                )
                cfile.write("\n\n")
                cfile.write(
                    "\n".join(
                        [
                            f"{tab}IJth(Jac, {idx//self.nspecies}, {idx%self.nspecies}) = {ccode(j)};"
                            for idx, j in enumerate(jac)
                        ]
                    )
                )
                # cfile.write(ccode(CodeBlock(*fex)))
                cfile.write("\n\n")
                cfile.write(f"{tab}return 0;")
                cfile.write("\n\n}")

    # TODO: combine with fex
    def jtv(self, symbol_form=False):
        rate_declare, rate_assign, rhs = self.rhs(symbol_form)

        ydot = MatrixSymbol("jv", self.nspecies, 1)
        fex_result = [Assignment(l, r) for l, r in zip(ydot, rhs)]
        return rate_declare, rate_assign, fex_result

    def jtv_to_ccode(
        self,
        *args,
        to_file=False,
        prefix="",
        fname="jtv",
        ext="c",
        header=False,
        h_ext="h",
        **kwargs,
    ):

        tab = "    "

        if to_file:

            cfunc = "int jtv(N_Vector v, N_Vector Jv, realtype t, N_Vector u, N_Vector fu, void *user_data, N_Vector tmp)"
            rate_declare, rate_assign, jtv = self.jtv(*args, **kwargs)

            with open(f"{prefix + fname}.{ext}", "w") as cfile:
                if header:
                    cfile.write(f'#include "{fname}.{h_ext}"\n\n')
                    with open(f"{prefix + fname}.{h_ext}", "w") as hfile:
                        hfile.write("#ifndef __JTV_H__\n")
                        hfile.write("#define __JTV_H__\n\n")
                        hfile.write("#include <cvode/cvode.h>\n")
                        hfile.write("#include <nvector/nvector_serial.h>\n\n")
                        hfile.write(f"{cfunc};\n\n")
                        hfile.write("#endif\n")
                else:
                    cfile.write("#include <cvode/cvode.h>\n")
                    cfile.write("#include <nvector/nvector_serial.h>\n\n")

                cfile.write("// Jacobian function vector routine.\n")
                cfile.write(f"{cfunc}")
                cfile.write("{\n\n")
                cfile.write(f"{tab}realtype *x  = N_VGetArrayPointer(u);\n")
                cfile.write(f"{tab}realtype *y  = N_VGetArrayPointer(v);\n")
                cfile.write(f"{tab}realtype *jv = N_VGetArrayPointer(Jv);\n")
                cfile.write(f"{tab}realtype *fx = N_VGetArrayPointer(fu);\n")
                cfile.write(f"{tab}UserData *u_data = (UserData*) user_data;\n")
                cfile.write(f"{tab}realtype Tgas = u_data->Tgas;\n\n")
                cfile.write(
                    "\n".join([f"{tab}{ccode(dec)} = 0.0;" for dec in rate_declare])
                )
                cfile.write("\n\n")
                cfile.write(
                    "\n".join(
                        [
                            f"{tab}if (Tgas>{react.temp_min} && Tgas<{react.temp_max}) {ccode(assign)}"
                            for react, assign in zip(self.reaction_list, rate_assign)
                        ]
                    )
                )
                cfile.write("\n\n")
                cfile.write("\n".join([tab + ccode(f) for f in jtv]))
                # cfile.write(ccode(CodeBlock(*fex)))
                cfile.write("\n\n")
                cfile.write(
                    "\n".join(f"{tab} fx[{i}] = 0.0;" for i in range(self.nspecies))
                )
                cfile.write("\n\n")
                cfile.write(f"{tab}return 0;")
                cfile.write("\n\n}")

        else:
            print("\n".join([ccode(f) for f in self.fex()]))

    def rhs(self, symbol_form=False):
        if not self.info_updated:
            self.collect_infos()

        y = MatrixSymbol("y", self.nspecies, 1)
        rhs_result = [0.0] * self.nspecies

        rate_declare = []
        rate_assign = []

        for r, react in enumerate(self.reaction_list):
            if symbol_form:
                # rate_symbol = symbols(f"rate_react{r}", cls=Function)()
                rate_symbol = Symbol(f"rate_react{r}")
                rate_declare.append(Declaration(Variable(rate_symbol, type=real)))
                rate_assign.append(Assignment(rate_symbol, react.rate_func()))
            else:
                rate_symbol = react.rate_func()
            ridx = [self.net_species.index(r) for r in react.reactants]
            pidx = [self.net_species.index(p) for p in react.products]
            rsym = [y[idx] for idx in ridx]
            rsym_product = reduce(mul, rsym)
            # psym = [y[idx] for idx in pidx]
            # psym_product = reduce(mul, psym)
            for idx in ridx:
                rhs_result[idx] -= rate_symbol * rsym_product
            for idx in pidx:
                rhs_result[idx] += rate_symbol * rsym_product

        return rate_declare, rate_assign, rhs_result
