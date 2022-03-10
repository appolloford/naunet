from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Type
from tqdm import tqdm
from jinja2 import Environment, PackageLoader

from .species import Species
from .reactions.reaction import Reaction
from .dusts.dust import Dust
from .utilities import _stmwrap

# class RelativeEnvironment(Environment):
#     """Override join_path() to enable relative template paths."""

#     def join_path(self, template, parent):
#         return os.path.join(parent, template)
#         # return os.path.join(os.path.dirname(parent), template)

# define in this file to avoid circular import
@dataclass
class NetworkInfo:
    n_spec: int
    n_react: int
    species: list[Species]
    reactions: list[Reaction]
    databases: list[Type[Reaction]]
    heating: list[str] = None
    cooling: list[str] = None
    dust: Dust = None
    shielding: dict = None


class TemplateLoader:
    """Class to render main templates

    The class receive the information provided by network and prepare
    the variables needed by the main templates (i.e. the required files
    to build the static/shared lib or python module of the chemical
    network).
    """

    @dataclass
    class InfoContent:
        """
        General solver information
        """

        method: str
        device: str

    @dataclass
    class MacrosContent:
        """
        The variables required by macros file
        """

        nspec: str
        neqns: str
        nreact: str
        nheating: str
        ncooling: str
        speclist: list[str]
        nnz: str = None

    @dataclass
    class ODEContent:
        """
        The information required by ode expressions
        """

        rateeqns: list[str]
        hrateeqns: list[str]
        crateeqns: list[str]
        fex: list[str]
        jac: list[str]
        jacpattern: list[str]
        spjacrptr: list[str] = None
        spjaccval: list[str] = None
        spjacrptrarr: str = None
        spjaccvalarr: str = None
        spjacdata: list[str] = None
        header: str = None
        odemodifier: list[str] = None
        ratemodifier: list[str] = None

    @dataclass
    class PhysicsContent:
        """
        The information required by general physics function
        """

        mantles: str
        Hnuclei: str
        mu: str
        gamma: str
        h2shielding: str = None
        coshielding: str = None
        n2shielding: str = None
        header: str = None

    @dataclass
    class VariablesContent:
        """
        The variables defined in reactions or thermal processes
        """

        consts: dict[str, str]
        globs: list[str]
        varis: list[str]
        user_var: list[str]
        tvaris: list[str]  # variables required by thermal process
        header: str = None

    def __init__(
        self,
        netinfo: NetworkInfo,
        solver: str,
        method: str,
        device: str,
        ratemodifier: dict[int, str] = None,
        odemodifier: list[str] = None,
    ) -> None:

        loader = PackageLoader("naunet")
        # self._env = RelativeEnvironment(loader=loader)
        self._env = Environment(loader=loader)
        self._env.globals.update(zip=zip)
        self._env.filters["stmwrap"] = _stmwrap
        self._solver = solver
        self._env.trim_blocks = True
        self._env.rstrip_blocks = True

        self._ode_modifier = odemodifier.copy() if odemodifier else []
        self._rate_modifier = ratemodifier.copy() if ratemodifier else {}

        self._info = None
        self._macros = None
        self._ode = None
        self._physics = None
        self._variables = None

        self._prepare_contents(netinfo, method, device)

    def _prepare_contents(self, netinfo: object, method: str, device: str) -> None:
        n_spec = netinfo.n_spec
        n_react = netinfo.n_react
        species = netinfo.species
        reactions = netinfo.reactions
        databases = netinfo.databases
        heating = netinfo.heating
        cooling = netinfo.cooling
        dust = netinfo.dust
        odemodifier = self._ode_modifier

        ratemodifier = [
            f"k[{idx}] = {value}"
            for idx, reac in enumerate(reactions)
            for key, value in self._rate_modifier.items()
            if key == reac.idxfromfile
        ]

        has_thermal = True if heating or cooling else False
        n_eqns = n_spec + has_thermal

        self._info = self.InfoContent(method, device)

        nspec = f"#define NSPECIES {n_spec}"
        neqns = f"#define NEQUATIONS {n_eqns}"
        nreact = f"#define NREACTIONS {n_react}"
        nheating = f"#define NHEATPROCS {len(heating) if heating else 0}"
        ncooling = f"#define NCOOLPROCS {len(cooling) if cooling else 0}"
        speclist = [f"#define IDX_{x.alias} {i}" for i, x in enumerate(species)]
        if has_thermal:
            speclist.append(f"#define IDX_TGAS {n_spec}")

        self._macros = self.MacrosContent(
            nspec,
            neqns,
            nreact,
            nheating,
            ncooling,
            speclist,
        )

        mantles = " + ".join(f"y[IDX_{g.alias}]" for g in species if g.is_surface)
        mantles = mantles if mantles else "0.0"

        Hnuclei = " + ".join(
            f"{s.element_count.get('H'):3.1e}*y[IDX_{s.alias}]"
            for s in species
            if "H" in s.element_count.keys()
        )
        Hnuclei = Hnuclei if Hnuclei else "0.0"

        # TODO: exclude electron, grain?
        density = " + ".join(f"y[IDX_{s.alias}]*{s.massnumber}" for s in species)
        npartile = " + ".join(f"y[IDX_{s.alias}]" for s in species)
        mu = "".join(["(", density, ") / (", npartile, ")"])

        # TODO: different ways to get adiabatic index
        gamma = "5.0 / 3.0"

        self._physics = self.PhysicsContent(
            mantles,
            Hnuclei,
            mu,
            gamma,
            h2shielding=netinfo.shielding.get("H2", ""),
            coshielding=netinfo.shielding.get("CO", ""),
            n2shielding=netinfo.shielding.get("N2", ""),
        )

        reactconsts = {
            f"{c:<15}": f"{cv}" for db in databases for c, cv in db.consts.items()
        }
        dustconsts = (
            {f"{c:<15}": f"{cv}" for c, cv in dust.consts.items()} if dust else {}
        )
        heatconsts = (
            {f"{c:<15}": f"{cv}" for p in heating for c, cv in p.consts.items()}
            if heating
            else {}
        )
        coolconsts = (
            {f"{c:<15}": f"{cv}" for p in cooling for c, cv in p.consts.items()}
            if cooling
            else {}
        )
        ebs = {
            f"eb_{s.alias:<12}": f"{s.binding_energy}" for s in species if s.is_surface
        }
        consts = {**reactconsts, **dustconsts, **heatconsts, **coolconsts, **ebs}

        reactglobs = [f"{v}" for db in databases for v in db.globs.values()]
        dustglobs = [f"{v}" for v in dust.globs.values()] if dust else []
        heatglobs = (
            [f"{v}" for p in heating for v in p.globs.values()] if heating else []
        )
        coolglobs = (
            [f"{v}" for p in cooling for v in p.globs.values()] if cooling else []
        )
        globs = [*reactglobs, *dustglobs, *heatglobs, *coolglobs]

        reactvars = [f"{v}" for db in databases for v in db.varis.values()]
        dustvars = [f"{v}" for v in dust.varis.values() if dust] if dust else []
        heatvars = (
            [f"{v}" for p in heating for v in p.varis.values()] if heating else []
        )
        coolvars = (
            [f"{v}" for p in cooling for v in p.varis.values()] if cooling else []
        )
        varis = [*dustvars, *reactvars]
        tvaris = ["mu", "gamma", *heatvars, *coolvars] if has_thermal else []

        react_uservar = [v for db in databases for v in db.user_var]
        dust_uservar = dust.user_var if dust else []
        heat_uservar = [v for p in heating for v in p.user_var]
        cool_uservar = [v for p in cooling for v in p.user_var]
        user_var = [*dust_uservar, *react_uservar, *heat_uservar, *cool_uservar]

        self._variables = self.VariablesContent(consts, globs, varis, user_var, tvaris)

        # prepare reaction rate expressions
        rates = [f"k[{r}]" for r in range(n_react)]
        rateeqns = [
            f"if (Tgas>{reac.temp_min} && Tgas<{reac.temp_max}) {{\n{' = '.join([rate, reac.rate_func()])}; \n}}"
            if reac.temp_min < reac.temp_max
            else f"{' = '.join([rate, reac.rate_func()])};"
            for rate, reac in zip(rates, reactions)
        ]

        y = [f"y[IDX_{x.alias}]" for x in species]
        if has_thermal:
            y.append("y[IDX_TGAS]")
        rhs = ["0.0"] * n_eqns
        jacrhs = ["0.0"] * n_eqns * n_eqns
        for rl, react in enumerate(tqdm(reactions, desc="Preparing ODE...")):

            rspecidx = [species.index(r) for r in react.reactants]
            pspecidx = [species.index(p) for p in react.products]

            # Differential Equation
            rsym = [y[idx] for idx in rspecidx]
            rsym_mul = "*".join(rsym)
            for specidx in rspecidx:
                rhs[specidx] += f" - {rates[rl]}*{rsym_mul}"
            for specidx in pspecidx:
                rhs[specidx] += f" + {rates[rl]}*{rsym_mul}"

            # Jacobian
            for specidx in rspecidx:
                # df/dx, remove the dependency for current reactant
                for ri in rspecidx:
                    rsymcopy = rsym.copy()
                    rsymcopy.remove(y[ri])
                    term = f" - {'*'.join([rates[rl], *rsymcopy])}"
                    jacrhs[specidx * n_eqns + ri] += term
            for specidx in pspecidx:
                for ri in rspecidx:
                    rsymcopy = rsym.copy()
                    rsymcopy.remove(y[ri])
                    term = f" + {'*'.join([rates[rl], *rsymcopy])}"
                    jacrhs[specidx * n_eqns + ri] += term

            # for specidx, _ in enumerate(species):

            #     # make the diagonal terms nonzero to make sure the matrix reversible
            #     if rhs[specidx] == "0.0":
            #         rhs[specidx] = "naunet_rate_tiny"

        # prepare heating rate expressions
        hrates = [f"kh[{hidx}]" for hidx, _ in enumerate(heating)]
        hrateeqns = [
            f"if (Tgas>{h.temp_min} && Tgas<{h.temp_max}) {{\n{' = '.join([hrate, h.rate_func()])}; \n}}"
            if h.temp_min < h.temp_max
            else f"{' = '.join([hrate, h.rate_func()])};"
            for hrate, h in zip(hrates, heating)
        ]

        # fex/jac of thermal process
        for hidx, h in enumerate(heating):

            rspecidx = [species.index(r) for r in h.reactants]
            rsym = [y[idx] for idx in rspecidx]
            rsym_mul = "*".join(rsym)

            # temperature index is n_spec
            rhs[n_spec] += f" + {hrates[hidx]} * {rsym_mul}"

            for ri in rspecidx:
                rsymcopy = rsym.copy()
                rsymcopy.remove(y[ri])
                term = f" + {'*'.join([hrates[hidx], *rsymcopy])}"
                # only fill the last row of jacobian
                jacrhs[n_spec * n_eqns + ri] += term

        # prepare cooling rate expressions
        crates = [f"kc[{cidx}]" for cidx, _ in enumerate(cooling)]
        crateeqns = [
            f"if (Tgas>{c.temp_min} && Tgas<{c.temp_max}) {{\n{' = '.join([crate, c.rate_func()])}; \n}}"
            if c.temp_min < c.temp_max
            else f"{' = '.join([crate, c.rate_func()])};"
            for crate, c in zip(crates, cooling)
        ]

        for cidx, c in enumerate(cooling):

            rspecidx = [species.index(r) for r in c.reactants]
            rsym = [y[idx] for idx in rspecidx]
            rsym_mul = "*".join(rsym)

            # temperature index is n_spec
            rhs[n_spec] += f" - {crates[cidx]} * {rsym_mul}"

            for ri in rspecidx:
                rsymcopy = rsym.copy()
                rsymcopy.remove(y[ri])
                term = f" - {'*'.join([crates[cidx], *rsymcopy])}"
                # only fill the last row of jacobian
                jacrhs[n_spec * n_eqns + ri] += term

        lhs = [f"ydot[IDX_{x.alias}]" for x in species]
        if has_thermal:
            lhs.append("ydot[IDX_TGAS]")
            rhs[n_spec] = f"(gamma - 1.0) * ( {rhs[n_spec]} ) / kerg / GetNumDens(y)"
            for si in range(n_spec):
                jacrhs[n_spec * n_eqns + si] = (
                    "0.0"
                    if jacrhs[n_spec * n_eqns + si] == "0.0"
                    else "".join(
                        [
                            "(gamma - 1.0) * (",
                            jacrhs[n_spec * n_eqns + si],
                            " ) / kerg / GetNumDens(y)",
                        ]
                    )
                )

        fex = [f"{l} = {r};" for l, r in zip(lhs, rhs)]

        if self._solver == "cvode":
            jac = [
                f"IJth(jmatrix, {idx//n_eqns}, {idx%n_eqns}) = {j};"
                for idx, j in enumerate(jacrhs)
                if j != "0.0"
            ]
        elif self._solver == "odeint":
            jac = [
                f"j({idx//n_eqns}, {idx%n_eqns}) = {j};"
                for idx, j in enumerate(jacrhs)
                if j != "0.0"
            ]

        jacpattern = [0 if j == "0.0" else 1 for j in jacrhs]

        rowpattern = []
        for row in range(n_eqns):
            rowdata = jacpattern[row * n_eqns : (row + 1) * n_eqns]
            rowpattern.append(" ".join(str(e) for e in rowdata))

        patternstr = "\n".join(rowpattern)

        spjacrptr = []
        spjaccval = []
        spjacdata = []
        spjacrptrarr = []
        spjaccvalarr = []
        if "sparse" in method:
            # TODO: Too slow! Optimize it
            nnz = 0
            for row in range(n_eqns):
                spjacrptr.append(f"rowptrs[{row}] = {nnz};")
                spjacrptrarr.append(str(nnz))
                for col in range(n_eqns):
                    elem = jacrhs[row * n_eqns + col]
                    if elem != "0.0":
                        spjaccval.append(f"colvals[{nnz}] = {col};")
                        spjaccvalarr.append(str(col))
                        spjacdata.append(f"data[{nnz}] = {elem};")
                        nnz += 1
            spjacrptr.append(f"rowptrs[{n_eqns}] = {nnz};")
            spjacrptrarr.append(str(nnz))
            self._macros.nnz = f"#define NNZ {nnz}"

        spjacrptrarr = ", ".join(spjacrptrarr)
        spjaccvalarr = ", ".join(spjaccvalarr)

        self._ode = self.ODEContent(
            rateeqns,
            hrateeqns,
            crateeqns,
            fex,
            jac,
            jacpattern=patternstr,
            spjacrptr=spjacrptr,
            spjaccval=spjaccval,
            spjacrptrarr=spjacrptrarr,
            spjaccvalarr=spjaccvalarr,
            spjacdata=spjacdata,
            odemodifier=odemodifier,
            ratemodifier=ratemodifier,
        )

    def _render(
        self, template: object, prefix: str, name: str, save: bool, *args, **kwargs
    ) -> None:

        # template = self.env.get_template(os.path.join(template_prefix, template_file))
        if save:
            template.stream(**kwargs).dump(os.path.join(prefix, name))

        else:
            result = template.render(**kwargs)
            print(result)

    def render(self, prefix="./source", save=True):

        Path(prefix).mkdir(parents=True, exist_ok=True)

        self.render_constants(prefix=prefix, save=save)
        self.render_macros(prefix=prefix, save=save)
        self.render_data(prefix=prefix, save=save)
        self.render_ode(prefix=prefix, save=save)
        self.render_physics(prefix=prefix, save=save)
        self.render_naunet(prefix=prefix, save=save)
        self.render_utilities(prefix=prefix, save=save)

    # This function should not be used outside console commands
    def render_cmake(self, prefix: str = "./", version: str = None) -> None:

        template_names = [
            "CMakeLists.txt.j2",
            "src/CMakeLists.txt.j2",
        ]
        for name in template_names:
            tname = os.path.join(self._solver, name)
            template = self._env.get_template(tname)
            target = name.replace(".j2", "")
            self._render(
                template,
                prefix,
                name=target,
                save=True,
                info=self._info,
                version=version,
            )

    def render_constants(
        self,
        prefix: str = "./",
        name: str = None,
        save: bool = True,
        headerprefix: str = None,
        headername: str = None,
        header: bool = True,
    ) -> None:

        if save:
            suffix = "cu" if self._info.device == "gpu" else "cpp"
            name = name if name else f"naunet_constants.{suffix}"
            headername = headername if headername else "naunet_constants.h"
            headerprefix = headerprefix if headerprefix else prefix

        if header:
            headername = headername if headername else "naunet_constants.h"
            self._variables.header = headername
            tname = "common/include/naunet_constants.h.j2"
            template = self._env.get_template(tname)
            self._render(
                template,
                headerprefix,
                headername,
                save,
                variables=self._variables,
                info=self._info,
                physics=self._physics,
            )

        tname = "common/src/naunet_constants.cpp.j2"
        template = self._env.get_template(tname)
        self._render(
            template,
            prefix,
            name,
            save,
            variables=self._variables,
            info=self._info,
            physics=self._physics,
        )

    def render_macros(
        self, prefix: str = "./", name: str = None, save: bool = True
    ) -> None:
        if not name and save:
            name = "naunet_macros.h"

        tname = "common/include/naunet_macros.h.j2"
        template = self._env.get_template(tname)
        self._render(template, prefix, name, save, macros=self._macros, info=self._info)

    def render_naunet(
        self,
        prefix: str = "./",
        name: str = None,
        save: bool = True,
        headerprefix: str = None,
        headername: str = None,
        header: bool = True,
    ):

        if save:
            name = name if name else "naunet.cpp"
            headername = headername if headername else "naunet.h"
            headerprefix = headerprefix if headerprefix else prefix

        if header:
            headername = headername if headername else "naunet.h"
            tname = os.path.join(self._solver, "include/naunet.h.j2")
            template = self._env.get_template(tname)
            self._render(
                template,
                headerprefix,
                headername,
                save,
                info=self._info,
                header=headername,
                variables=self._variables,
            )

        tname = os.path.join(self._solver, "src/naunet.cpp.j2")
        template = self._env.get_template(tname)
        self._render(
            template,
            prefix,
            name,
            save,
            info=self._info,
            header=headername,
            physics=self._physics,
            variables=self._variables,
        )

    def render_ode(
        self,
        prefix: str = "./",
        name: str = None,
        save: bool = True,
        headerprefix: str = None,
        headername: str = None,
        header: bool = True,
    ) -> None:

        if save:
            suffix = "cu" if self._info.device == "gpu" else "cpp"
            name = name if name else f"naunet_ode.{suffix}"
            headername = headername if headername else "naunet_ode.h"
            headerprefix = headerprefix if headerprefix else prefix

        if header:
            headername = headername if headername else "naunet_ode.h"
            self._ode.header = headername
            tname = os.path.join(self._solver, "include/naunet_ode.h.j2")
            template = self._env.get_template(tname)
            self._render(
                template,
                headerprefix,
                headername,
                save,
                ode=self._ode,
                info=self._info,
            )

        if self._solver == "cvode":
            tname = os.path.join(self._solver, "src/naunet_rates.cpp.j2")
            template = self._env.get_template(tname)
            self._render(
                template,
                prefix,
                f"naunet_rates.{suffix}",
                save,
                ode=self._ode,
                info=self._info,
                variables=self._variables,
            )
            tname = os.path.join(self._solver, "src/naunet_fex.cpp.j2")
            template = self._env.get_template(tname)
            self._render(
                template,
                prefix,
                f"naunet_fex.{suffix}",
                save,
                ode=self._ode,
                info=self._info,
                variables=self._variables,
            )
            tname = os.path.join(self._solver, "src/naunet_jac.cpp.j2")
            template = self._env.get_template(tname)
            self._render(
                template,
                prefix,
                f"naunet_jac.{suffix}",
                save,
                ode=self._ode,
                info=self._info,
                variables=self._variables,
            )

        else:
            tname = os.path.join(self._solver, "src/naunet_ode.cpp.j2")
            template = self._env.get_template(tname)
            self._render(
                template,
                prefix,
                name,
                save,
                ode=self._ode,
                info=self._info,
                variables=self._variables,
            )

    def render_physics(
        self,
        prefix: str = "./",
        name: str = None,
        save: bool = True,
        headerprefix: str = None,
        headername: str = None,
        header: bool = True,
    ):

        if save:
            suffix = "cu" if self._info.device == "gpu" else "cpp"
            name = name if name else f"naunet_physics.{suffix}"
            headername = headername if headername else "naunet_physics.h"
            headerprefix = headerprefix if headerprefix else prefix

        if header:
            headername = headername if headername else "naunet_physics.h"
            self._physics.header = headername
            tname = "common/include/naunet_physics.h.j2"
            template = self._env.get_template(tname)
            self._render(
                template,
                headerprefix,
                headername,
                save,
                physics=self._physics,
                info=self._info,
            )

        tname = "common/src/naunet_physics.cpp.j2"
        template = self._env.get_template(tname)
        self._render(
            template,
            prefix,
            name,
            save,
            physics=self._physics,
            info=self._info,
        )

    def render_data(
        self, prefix: str = "./", name: str = None, save: bool = True
    ) -> None:
        if not name and save:
            name = "naunet_data.h"

        tname = "common/include/naunet_data.h.j2"
        template = self._env.get_template(tname)
        self._render(template, prefix, name, save, variables=self._variables)

    def render_jac_pattern(self, prefix: str = "./") -> None:

        with open(f"{prefix}/jac_pattern.dat", "w") as outf:
            outf.write(self._ode.jacpattern)

    def render_utilities(
        self,
        prefix: str = "./",
        name: str = None,
        save: bool = True,
        headerprefix: str = None,
        headername: str = None,
        header: bool = True,
    ):

        if save:
            name = name if name else f"naunet_utilities.cpp"
            headername = headername if headername else "naunet_utilities.h"
            headerprefix = headerprefix if headerprefix else prefix

        if header:
            headername = headername if headername else "naunet_utilities.h"
            self._physics.header = headername
            tname = "common/include/naunet_utilities.h.j2"
            template = self._env.get_template(tname)
            self._render(template, headerprefix, headername, save)

        tname = "common/src/naunet_utilities.cpp.j2"
        template = self._env.get_template(tname)
        self._render(template, prefix, name, save)
