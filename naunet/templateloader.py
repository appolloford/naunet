import os
from dataclasses import dataclass
from pathlib import Path
from typing import List, Dict
from tqdm import tqdm
from jinja2 import Environment, PackageLoader


# class RelativeEnvironment(Environment):
#     """Override join_path() to enable relative template paths."""

#     def join_path(self, template, parent):
#         return os.path.join(parent, template)
#         # return os.path.join(os.path.dirname(parent), template)


class TemplateLoader:
    @dataclass
    class InfoContent:
        method: str
        device: str

    @dataclass
    class MacrosContent:
        nspec: str
        nreact: str
        speclist: List[str]
        nnz: str = None

    @dataclass
    class ODEContent:
        method: str
        device: str
        rateeqns: List[str]
        fex: List[str]
        jac: List[str]
        spjacrptr: List[str] = None
        spjaccval: List[str] = None
        spjacdata: List[str] = None
        header: str = None
        odemodifier: List[str] = None
        ratemodifier: List[str] = None

    @dataclass
    class PhysicsContent:
        mantles: str
        h2shielding: str = None
        coshielding: str = None
        n2shielding: str = None
        header: str = None

    @dataclass
    class VariablesContent:
        consts: Dict[str, str]
        globs: List[str]
        vars: List[str]
        user_var: List[str]
        header: str = None

    def __init__(self, netinfo: object, solver: str, method: str, device: str) -> None:

        loader = PackageLoader("naunet")
        # self._env = RelativeEnvironment(loader=loader)
        self._env = Environment(loader=loader)
        self._env.globals.update(zip=zip)
        self._solver = solver
        self._env.trim_blocks = True
        self._env.rstrip_blocks = True

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
        dust = netinfo.dust
        odemodifier = netinfo.odemodifier
        ratemodifier = netinfo.ratemodifier

        self._info = self.InfoContent(method, device)

        nspec = f"#define NSPECIES {n_spec}"
        nreact = f"#define NREACTIONS {n_react}"
        speclist = [f"#define IDX_{x.alias} {i}" for i, x in enumerate(species)]

        self._macros = self.MacrosContent(nspec, nreact, speclist)

        mantles = " + ".join(f"y[IDX_{g.alias}]" for g in species if g.is_surface)
        mantles = mantles if mantles else "0.0"
        self._physics = self.PhysicsContent(
            mantles,
            h2shielding="L96Table",
            coshielding="V09Table",
            n2shielding="L13Table",
        )

        reactconsts = {
            f"{c:<15}": f"{cv}" for db in databases for c, cv in db.consts.items()
        }
        dustconsts = (
            {f"{c:<15}": f"{cv}" for c, cv in dust.consts.items()} if dust else {}
        )
        ebs = {
            f"eb_{s.alias:<12}": f"{s.binding_energy}" for s in species if s.is_surface
        }
        consts = {**reactconsts, **dustconsts, **ebs}

        reactglobs = [f"{v}" for db in databases for v in db.vars.values()]
        dustglobs = [f"{v}" for v in dust.globs.values()] if dust else []
        globs = [*reactglobs, *dustglobs]

        reactvars = [f"{v}" for db in databases for v in db.vars.values()]
        dustvars = [f"{v}" for v in dust.vars.values() if dust] if dust else []
        vars = [*dustvars, *reactvars]

        react_uservar = [v for db in databases for v in db.user_var]
        dust_uservar = dust.user_var if dust else []
        user_var = [*dust_uservar, *react_uservar]

        self._variables = self.VariablesContent(consts, globs, vars, user_var)

        rates = [f"k[{r}]" for r in range(n_react)]
        rateeqns = [
            f"if (Tgas>{reac.temp_min} && Tgas<{reac.temp_max}) {{\n{' = '.join([rate, reac.rate_func()])}; \n}}"
            if reac.temp_min < reac.temp_max
            else f"{' = '.join([rate, reac.rate_func()])};"
            for rate, reac in zip(rates, reactions)
        ]

        y = [f"y[IDX_{x.alias}]" for x in species]
        rhs = ["0.0"] * n_spec
        jacrhs = ["0.0"] * n_spec * n_spec
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
                    jacrhs[specidx * n_spec + ri] += term
            for specidx in pspecidx:
                for ri in rspecidx:
                    rsymcopy = rsym.copy()
                    rsymcopy.remove(y[ri])
                    term = f" + {'*'.join([rates[rl], *rsymcopy])}"
                    jacrhs[specidx * n_spec + ri] += term

        lhs = [f"ydot[IDX_{x.alias}]" for x in species]
        fex = [f"{l} = {r};" for l, r in zip(lhs, rhs)]

        if self._solver == "cvode":
            jac = [
                f"IJth(jmatrix, {idx//n_spec}, {idx%n_spec}) = {j};"
                for idx, j in enumerate(jacrhs)
            ]
        elif self._solver == "odeint":
            jac = [
                f"j({idx//n_spec}, {idx%n_spec}) = {j};" for idx, j in enumerate(jacrhs)
            ]

        spjacrptr = []
        spjaccval = []
        spjacdata = []
        if "sparse" in method:
            # TODO: Too slow! Optimize it
            nnz = 0
            for row in range(n_spec):
                spjacrptr.append(f"rowptrs[{row}] = {nnz};")
                for col in range(n_spec):
                    elem = jacrhs[row * n_spec + col]
                    if elem != "0.0":
                        spjaccval.append(f"colvals[{nnz}] = {col};")
                        spjacdata.append(f"data[{nnz}] = {elem};")
                        nnz += 1
            spjacrptr.append(f"rowptrs[{n_spec}] = {nnz};")
            self._macros.nnz = f"#define NNZ {nnz}"

        self._ode = self.ODEContent(
            method,
            device,
            rateeqns,
            fex,
            jac,
            spjacrptr=spjacrptr,
            spjaccval=spjaccval,
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
        self.render_userdata(prefix=prefix, save=save)
        self.render_ode(prefix=prefix, save=save)
        self.render_physics(prefix=prefix, save=save)
        self.render_naunet(prefix=prefix, save=save)

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
            name = (
                name
                if name
                else "naunet_constants.cu"
                if self._info.device == "gpu"
                else "naunet_constants.cpp"
            )
            headername = headername if headername else "naunet_constants.h"
            headerprefix = headerprefix if headerprefix else prefix

        if header:
            headername = headername if headername else "naunet_constants.h"
            self._variables.header = headername
            tname = os.path.join(self._solver, "include/naunet_constants.h.j2")
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

        tname = os.path.join(self._solver, "src/naunet_constants.cpp.j2")
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

        tname = os.path.join(self._solver, "include/naunet_macros.h.j2")
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
            )

        tname = os.path.join(self._solver, "src/naunet.cpp.j2")
        template = self._env.get_template(tname)
        self._render(template, prefix, name, save, info=self._info, header=headername)

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
            name = (
                name
                if name
                else "naunet_ode.cu"
                if self._info.device == "gpu"
                else "naunet_ode.cpp"
            )
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
            name = (
                name
                if name
                else "naunet_physics.cu"
                if self._info.device == "gpu"
                else "naunet_physics.cpp"
            )
            headername = headername if headername else "naunet_physics.h"
            headerprefix = headerprefix if headerprefix else prefix

        if header:
            headername = headername if headername else "naunet_physics.h"
            self._physics.header = headername
            tname = os.path.join(self._solver, "include/naunet_physics.h.j2")
            template = self._env.get_template(tname)
            self._render(
                template,
                headerprefix,
                headername,
                save,
                physics=self._physics,
                info=self._info,
            )

        tname = os.path.join(self._solver, "src/naunet_physics.cpp.j2")
        template = self._env.get_template(tname)
        self._render(
            template,
            prefix,
            name,
            save,
            physics=self._physics,
            info=self._info,
        )

    def render_userdata(
        self, prefix: str = "./", name: str = None, save: bool = True
    ) -> None:
        if not name and save:
            name = "naunet_userdata.h"

        tname = os.path.join(self._solver, "include/naunet_userdata.h.j2")
        template = self._env.get_template(tname)
        self._render(template, prefix, name, save, variables=self._variables)
