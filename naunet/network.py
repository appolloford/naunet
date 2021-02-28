import logging
from operator import mul
from functools import reduce
from sympy import symbols, Symbol, Function, MatrixSymbol, Idx, ccode
from sympy.codegen.ast import Assignment, CodeBlock, Declaration, Variable, real
from sympy.utilities.codegen import codegen
from jinja2 import Environment, FileSystemLoader, PackageLoader
from .settings import ode_symbols, user_symbols
from .species import Species
from .reactions.reaction import Reaction
from .reactions.kidareaction import KIDAReaction
from .reactions.leedsreaction import LEEDSReaction
from .reactions.kromereaction import KROMEReaction

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()


supported_reaction_class = {
    "kida": KIDAReaction,
    "leeds": LEEDSReaction,
    "krome": KROMEReaction,
}


krome_globals = {
    "format": "idx,r,r,r,p,p,p,p,tmin,tmax,rate",
}


def krome_parser(line):
    if line.startswith(("#", "//")):
        return ""
    elif line.startswith("@format:"):
        krome_globals["format"] = line.replace("@format:", "")
        return ""
    else:
        return line.strip()


def reaction_factory(react_string: str, database: str) -> Reaction:
    initializer = supported_reaction_class.get(database)
    format = krome_globals.get("format")
    return initializer(react_string, format)


class TemplateLoader:
    def __init__(self) -> None:
        rpath = "cxx_src/cvode_example"
        self.loader = PackageLoader("naunet", rpath)
        self.env = Environment(loader=self.loader)
        self.env.trim_blocks = True
        # env.lstrip_blocks = True
        self.env.rstrip_blocks = True

    def _write_code(
        self,
        template_prefix: str,
        template_file: str,
        kw_write: list,
        *args,
        use_template: bool = True,
        to_file: bool = True,
        file_name: str = None,
        prefix: str = None,
        **kwargs,
    ):

        if use_template:
            template = self.env.get_template(template_prefix + template_file)

            if to_file:
                template.stream(**kwargs).dump(prefix + file_name)

            else:
                result = template.render()
                print(result)

        else:

            if len(kw_write):

                if to_file:
                    with open(prefix + file_name, "w") as hfile:
                        for kw in kw_write:
                            value = kwargs.get(kw)
                            if isinstance(value, list):
                                outstr = "\n".join(str(x) for x in value)
                            else:
                                outstr = str(value)

                            hfile.write(f"{outstr}\n")

                else:
                    for kw in kw_write:
                        value = kwargs.get(kw)
                        if isinstance(value, list):
                            outstr = "\n".join(str(x) for x in value)
                        else:
                            outstr = str(value)
                        print(f"{outstr}\n")


class Info(TemplateLoader):
    def __init__(self, speclist: list) -> None:

        super().__init__()

        self.net_species = speclist
        self.n_spec = len(speclist)

    # ? create the to_code() interface in TemplateLoader
    def to_ccode(self, *args, **kwargs) -> None:

        nspec_def = f"#define NSPECIES {self.n_spec}\n"

        spec_idx = [
            f"#define IDX_{x.alias} {i}" for i, x in enumerate(self.net_species)
        ]

        template_prefix = "include/"
        template_file = "naunet_constants.h.j2"
        kw_write = ["nspec_def", "spec_idx"]

        self._write_code(
            template_prefix,
            template_file,
            kw_write,
            **kwargs,
            nspec_def=nspec_def,
            spec_idx=spec_idx,
        )

        # if use_template:
        #     template = self.env.get_template(f"include/naunet_constants.h.j2")

        #     if to_file:
        #         template.stream(nspec_def=nspec_def, spec_idx=spec_idx).dump(
        #             f"{prefix}{filename}"
        #         )

        #     else:
        #         result = template.render()
        #         print(result)

        # else:

        #     if to_file:
        #         with open(f"{prefix}{filename}", "w") as hfile:
        #             hfile.write(f"#define NSPECIES {self.n_spec}\n\n")
        #             hfile.write("\n".join(spec_idx))

        #     else:
        #         print("\n".join(spec_idx))


class UserData(TemplateLoader):
    def __init__(self) -> None:
        super().__init__()

    def to_ccode(self, *args, **kwargs) -> None:

        variables = [
            ccode(Declaration(Variable(v, type=real))) for v in user_symbols.values()
        ]

        template_prefix = "include/"
        template_file = "naunet_userdata.h.j2"
        kw_write = ["variables"]

        self._write_code(
            template_prefix, template_file, kw_write, **kwargs, variables=variables
        )


class Network:
    def __init__(self, fname: str = None, database: str = None) -> None:
        self.reaction_list = []
        self.reactants_in_network = set()
        self.products_in_network = set()
        self.info = None
        self.userdata = UserData()
        self.ode_expression = None

        self.add_reaction_from_file(fname, database)

    def _add_reaction(self, react_string: str, database: str) -> list:
        reaction = reaction_factory(react_string, database)
        self.reaction_list.append(reaction)
        new_reactants = set(reaction.reactants).difference(self.reactants_in_network)
        new_products = set(reaction.products).difference(self.products_in_network)
        self.reactants_in_network.update(new_reactants)
        self.products_in_network.update(new_products)
        if len(self.reaction_list) % 100 == 0:
            print("Processing: {} reactions...".format(len(self.reaction_list)))
        return new_reactants | new_products

    def add_reaction(self, react_string: str, database: str) -> None:
        new_species = self._add_reaction(react_string, database)
        logger.info("New species are added: {}".format(new_species))

        # reset network information if content is changed
        self.info = None

    def add_reaction_from_file(self, filename: str, database: str) -> None:
        if not database:
            logger.critical(
                'Try to read in file but database is not assigned. Try again by "add_reaction_from_file"'
            )
        new_species = set()
        with open(filename, "r") as networkfile:
            for line in networkfile.readlines():
                if database == "krome":
                    react_string = krome_parser(line)
                    if react_string:
                        new_species.update(self._add_reaction(line, "krome"))
                else:
                    new_species.update(self._add_reaction(line, database))

            print("New species: \n{}".format("\n".join(str(x) for x in new_species)))

        self.info = None

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

    def get_info(self):
        if self.info:
            return self.info

        speclist = list(self.reactants_in_network | self.products_in_network)
        self.info = Info(speclist)
        # self.net_species = list(self.reactants_in_network | self.products_in_network)
        # self.nspecies = len(self.net_species)
        logging.info(
            "{} species in the network: {}".format(
                self.info.n_spec, ", ".join([x.name for x in self.info.net_species])
            )
        )

        # if information is re-collected, ode system must be reset
        self.ode_expression = None
        return self.info

    def create_ode_expression(self):

        # return the saved ode system if it has been updated
        if self.info and self.ode_expression:
            return self.ode_expression

        if not self.info:
            self.get_info()

        # renew an ode system
        self.ode_expression = ODESystem(self.info.n_spec, len(self.reaction_list))

        y = self.ode_expression.y
        rate_sym = self.ode_expression.rate_sym

        for rl, react in enumerate(self.reaction_list):
            # self.ode_expression.rate_func[rl] = react.rate_func()
            self.ode_expression.rate_func.append(react.rate_func())
            self.ode_expression.rate_mintemp.append(react.temp_min)
            self.ode_expression.rate_maxtemp.append(react.temp_max)

            ridx = [self.info.net_species.index(r) for r in react.reactants]
            pidx = [self.info.net_species.index(p) for p in react.products]

            rsym = [y[idx] for idx in ridx]
            rsym_mul = reduce(mul, rsym)
            for idx in ridx:
                self.ode_expression.rhs[idx] -= rate_sym[rl] * rsym_mul
            for idx in pidx:
                self.ode_expression.rhs[idx] += rate_sym[rl] * rsym_mul

            for idx in ridx:
                rsym_mul = reduce(mul, rsym)
                for ri in ridx:
                    residue = rsym_mul / y[ri]
                    self.ode_expression.jac[idx * self.info.n_spec + ri] -= (
                        rate_sym[rl] * residue
                    )
            for idx in pidx:
                rsym_mul = reduce(mul, rsym)
                for ri in ridx:
                    residue = rsym_mul / y[ri]
                    self.ode_expression.jac[idx * self.info.n_spec + ri] += (
                        rate_sym[rl] * residue
                    )

        return self.ode_expression


class ODESystem(TemplateLoader):

    cpp_func_names = {
        "cvode_fex": "int fex(realtype t, N_Vector u, N_Vector u_dot, void *user_data)",
        "cvode_jtv": "int jtv(N_Vector v, N_Vector Jv, realtype t, N_Vector u, N_Vector fu, void *user_data, N_Vector tmp)",
        "cvode_jac": "int jac(realtype t, N_Vector u, N_Vector fu, SUNMatrix Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)",
    }

    def __init__(self, n_spec: int, n_react: int) -> None:

        super().__init__()

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
        *args,
        var_rate: bool = True,
        header: bool = True,
        header_file: str = None,
        solver: str = "cvode",
        device: str = "cpu",
        **kwargs,
    ):

        cpp_func = self.cpp_func_names[f"{solver}_{function}"]

        rate_declare = [
            ccode(Declaration(Variable(s, type=real))) for s in self.rate_sym
        ]

        rate_assign = [
            f"if (Tgas>{tmin} && Tgas<{tmax}) {{\n{ccode(Assignment(sym, func))} \n}}"
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

        template_prefix = "src/"
        template_file = f"cv_{function}.cpp.j2"
        kw_write = ["rate_declare", "rate_assign", "eqns"]

        self._write_code(
            template_prefix,
            template_file,
            kw_write,
            **kwargs,
            header=header,
            header_file=header_file,
            func=cpp_func,
            vector=ode_symbols["ode_vector"],
            lhs=ode_symbols[f"{function}_lhs"],
            rate_declare=rate_declare,
            rate_assign=rate_assign,
            eqns=eqns,
        )

        if header:
            if not header_file:
                logger.warning(
                    'Header is used but file name is not assigned, "use header_file="'
                )

            template_prefix = "include/"
            template_file = f"cv_{function}.h.j2"
            kw_write = []

            self._write_code(
                template_prefix,
                template_file,
                kw_write,
                to_file=kwargs.get("to_file"),
                file_name=header_file,
                prefix=kwargs.get("prefix"),
                header=header,
                func=cpp_func,
            )

        # if use_template:

        #     template = self.env.get_template(f"src/cv_{function}.cpp.j2")
        #     if to_file:
        #         template.stream(
        #             header=header,
        #             file_name=file_name,
        #             h_ext=h_ext,
        #             func=cpp_func,
        #             vector=ode_symbols["ode_vector"],
        #             lhs=ode_symbols[f"{function}_lhs"],
        #             rate_declare=rate_declare,
        #             rate_assign=rate_assign,
        #             eqns=eqns,
        #         ).dump(f"{prefix + file_name}.{ext}")
        #     else:
        #         output = template.render(
        #             header=header,
        #             func=cpp_func,
        #             vector=ode_symbols["ode_vector"],
        #             lhs=ode_symbols[f"{function}_lhs"],
        #             rate_declare=rate_declare,
        #             rate_assign=rate_assign,
        #             eqns=eqns,
        #         )
        #         print(output)

        #     if header:
        #         template = self.env.get_template(f"include/cv_{function}.h.j2")
        #         if to_file:
        #             template.stream(func=cpp_func).dump(f"{prefix + file_name}.{h_ext}")
        #         else:
        #             output = template.render(fex_func_name=cpp_func)
        #             print(output)
        # else:

        #     if to_file:
        #         with open(f"{prefix + file_name}.{ext}", "w") as cfile:
        #             cfile.write("\n".join(ccode(dec) for dec in rate_declare))
        #             cfile.write("\n".join(ccode(assign) for assign in rate_declare))
        #             cfile.write("\n".join(ccode(eq) for eq in eqns))
        #     else:
        #         print(ccode(dec) for dec in rate_declare)
        #         print(ccode(assign) for assign in rate_declare)
        #         print("\n".join(ccode(eq) for eq in eqns))
