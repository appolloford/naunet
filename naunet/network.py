import logging
from operator import mul
from functools import reduce
from sympy import symbols, Symbol, Function, MatrixSymbol, Idx, ccode
from sympy.codegen.ast import Assignment, CodeBlock, Declaration, Variable, real
from sympy.utilities.codegen import codegen
from jinja2 import Environment, FileSystemLoader, PackageLoader
from .settings import ode_symbols
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
        self.ode_expression = None

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
        # if information is re-collected, ode system must be reset
        self.ode_expression = None

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

    # def ode_expr_func(self, form: str, var_rate: bool = True):
    #     if not self.info_updated:
    #         self.collect_infos()

    #     y = MatrixSymbol("y", self.nspecies, 1)
    #     if form == "rhs":
    #         expr_func = [0.0] * self.nspecies
    #     elif form == "jac":
    #         expr_func = [0.0] * self.nspecies * self.nspecies
    #     else:
    #         raise ValueError(f"Unknown value of form: {form}")

    #     rate_declare = []
    #     rate_assign = []

    #     for r, react in enumerate(self.reaction_list):
    #         if var_rate:
    #             # rate_symbol = symbols(f"rate_react{r}", cls=Function)()
    #             rate_symbol = Symbol(f"rate_react{r}")
    #             rate_declare.append(Declaration(Variable(rate_symbol, type=real)))
    #             rate_assign.append(Assignment(rate_symbol, react.rate_func()))
    #         else:
    #             rate_symbol = react.rate_func()

    #         ridx = [self.net_species.index(r) for r in react.reactants]
    #         pidx = [self.net_species.index(p) for p in react.products]

    #         if form == "rhs":

    #             rsym = [y[idx] for idx in ridx]
    #             rsym_product = reduce(mul, rsym)
    #             # psym = [y[idx] for idx in pidx]
    #             # psym_product = reduce(mul, psym)
    #             for idx in ridx:
    #                 expr_func[idx] -= rate_symbol * rsym_product
    #             for idx in pidx:
    #                 expr_func[idx] += rate_symbol * rsym_product

    #         elif form == "jac":

    #             for idx in ridx:
    #                 for ri in ridx:
    #                     rsym = [y[i] for i in ridx if i != ri]
    #                     rsym_product = reduce(mul, rsym)
    #                     expr_func[idx * self.nspecies + ri] -= (
    #                         rate_symbol * rsym_product
    #                     )
    #             for idx in pidx:
    #                 for ri in ridx:
    #                     rsym = [y[i] for i in ridx if i != ri]
    #                     rsym_product = reduce(mul, rsym)
    #                     expr_func[idx * self.nspecies + ri] += (
    #                         rate_symbol * rsym_product
    #                     )

    #     return rate_declare, rate_assign, expr_func

    # def ode_expression(self, func_name: str, *args, **kwargs):

    #     # TODO: Return a pure sympy object list when it is possible, make it can be handle by ccode()
    #     # * Sympy doesn't generate semicolons after declaration, neither support the if statement.
    #     # * Therefore, ode_expression() return the declaration and assignment with the final function.
    #     # * expr_to_code() will do further processing to generate the correct code.
    #     # fex_result = rate_declare
    #     # fex_result += rate_assign

    #     rate_declare, rate_assign, ode_expr_func = self.ode_expr_func(*args, **kwargs)

    #     if func_name == "fex":
    #         lhs = MatrixSymbol("ydot", self.nspecies, 1)
    #         expression = [Assignment(l, r) for l, r in zip(lhs, ode_expr_func)]
    #     elif func_name == "jtv":
    #         lhs = MatrixSymbol("jv", self.nspecies, 1)
    #         expression = [Assignment(l, r) for l, r in zip(lhs, ode_expr_func)]
    #     elif func_name == "jac":
    #         expression = ode_expr_func
    #     else:
    #         raise ValueError(f"Unknown value func_name: {func_name}")

    #     return rate_declare, rate_assign, expression

    def create_ode_expression(self):

        # return the saved ode system if it has been updated
        if self.info_updated and self.ode_expression:
            return self.ode_expression

        if not self.info_updated:
            self.collect_infos()

        # renew an ode system
        self.ode_expression = ODESystem(self.nspecies, len(self.reaction_list))

        y = self.ode_expression.y
        rate_sym = self.ode_expression.rate_sym

        for rl, react in enumerate(self.reaction_list):
            # self.ode_expression.rate_func[rl] = react.rate_func()
            self.ode_expression.rate_func.append(react.rate_func())
            self.ode_expression.rate_mintemp.append(react.temp_min)
            self.ode_expression.rate_maxtemp.append(react.temp_max)

            ridx = [self.net_species.index(r) for r in react.reactants]
            pidx = [self.net_species.index(p) for p in react.products]

            rsym = [y[idx] for idx in ridx]
            rsym_mul = reduce(mul, rsym)
            for idx in ridx:
                self.ode_expression.rhs[idx] -= rate_sym[rl] * rsym_mul
            for idx in pidx:
                self.ode_expression.rhs[idx] += rate_sym[rl] * rsym_mul

            for idx in ridx:
                for ri in ridx:
                    rsym = [y[i] for i in ridx if i != ri]
                    rsym_mul = reduce(mul, rsym)
                    self.ode_expression.jac[idx * self.nspecies + ri] -= (
                        rate_sym[rl] * rsym_mul
                    )
            for idx in pidx:
                for ri in ridx:
                    rsym = [y[i] for i in ridx if i != ri]
                    rsym_mul = reduce(mul, rsym)
                    self.ode_expression.jac[idx * self.nspecies + ri] += (
                        rate_sym[rl] * rsym_mul
                    )

        return self.ode_expression


class ODESystem:

    file_loader = PackageLoader("naunet", "cxx_src/cvode_example")
    env = Environment(loader=file_loader)
    env.trim_blocks = True
    # env.lstrip_blocks = True
    env.rstrip_blocks = True

    cpp_func_collect = {
        "cvode_fex": "int fex(realtype t, N_Vector u, N_Vector u_dot, void *user_data)",
        "cvode_jtv": "int jtv(N_Vector v, N_Vector Jv, realtype t, N_Vector u, N_Vector fu, void *user_data, N_Vector tmp)",
        "cvode_jac": "int jac(realtype t, N_Vector u, N_Vector fu, SUNMatrix Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)",
    }

    def __init__(self, n_spec: int, n_react: int) -> None:
        self.neq = n_spec
        self.nreact = n_react

        self.rate_sym = [Symbol(f"{ode_symbols['rate']}{r}") for r in range(n_react)]
        self.rate_func = []
        self.rate_mintemp = []
        self.rate_maxtemp = []

        self.y = MatrixSymbol(ode_symbols["ode_vector"], n_spec, 1)
        self.rhs = [0.0] * n_spec
        self.jac = [0.0] * n_spec * n_spec

    def to_ccode(
        self,
        function: str,
        var_rate: bool = True,
        use_template: bool = True,
        to_file: bool = True,
        file_name: str = None,
        prefix: str = None,
        ext: str = "cpp",
        header: bool = True,
        h_ext: str = "h",
        solver: str = "cvode",
        device: str = "cpu",
    ):

        cpp_func = self.cpp_func_collect[f"{solver}_{function}"]

        rate_declare = [
            ccode(Declaration(Variable(s, type=real))) for s in self.rate_sym
        ]

        rate_assign = [
            f"if (Tgas>{tmin} && Tgas<{tmax}) {ccode(Assignment(sym, func))}"
            if tmin < tmax
            else ccode(Assignment(sym, func))
            for tmin, tmax, sym, func in zip(
                self.rate_mintemp, self.rate_maxtemp, self.rate_sym, self.rate_func
            )
        ]

        if function == "jac":
            eqns = [
                f"IJth(Jac, {idx//self.neq}, {idx%self.neq}) = {ccode(j)};"
                for idx, j in enumerate(self.jac)
            ]
        else:
            lhs = MatrixSymbol(ode_symbols[f"{function}_lhs"], self.neq, 1)
            eqns = [ccode(Assignment(l, r)) for l, r in zip(lhs, self.rhs)]

        if use_template:

            template = self.env.get_template(f"src/cv_{function}.cpp.j2")
            if to_file:
                template.stream(
                    header=header,
                    file_name=file_name,
                    h_ext=h_ext,
                    func=cpp_func,
                    vector=ode_symbols["ode_vector"],
                    lhs=ode_symbols[f"{function}_lhs"],
                    rate_declare=rate_declare,
                    rate_assign=rate_assign,
                    eqns=eqns,
                ).dump(f"{prefix + file_name}.{ext}")
            else:
                output = template.render(
                    header=header,
                    func=cpp_func,
                    vector=ode_symbols["ode_vector"],
                    lhs=ode_symbols[f"{function}_lhs"],
                    rate_declare=rate_declare,
                    rate_assign=rate_assign,
                    eqns=eqns,
                )
                print(output)

            if header:
                template = self.env.get_template(f"include/cv_{function}.h.j2")
                if to_file:
                    template.stream(func=cpp_func).dump(f"{prefix + file_name}.{h_ext}")
                else:
                    output = template.render(fex_func_name=cpp_func)
                    print(output)
        else:

            if to_file:
                with open(f"{prefix + file_name}.{ext}", "w") as cfile:
                    cfile.write("\n".join(ccode(dec) for dec in rate_declare))
                    cfile.write("\n".join(ccode(assign) for assign in rate_declare))
                    cfile.write("\n".join(ccode(eq) for eq in eqns))
            else:
                print(ccode(dec) for dec in rate_declare)
                print(ccode(assign) for assign in rate_declare)
                print("\n".join(ccode(eq) for eq in eqns))
