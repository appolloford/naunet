import os
from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import List
from tqdm import tqdm
from jinja2 import Environment, FileSystemLoader, PackageLoader
from . import settings

# @dataclass
# class AbstractDataclass(ABC):
#     def __new__(cls, *args, **kwargs):
#         if cls == AbstractDataclass or cls.__bases__[0] == AbstractDataclass:
#             raise TypeError("Cannot instantiate abstract class.")
#         return super().__new__(cls)


class TemplateLoader(ABC):
    @dataclass
    class MacrosContent:
        nspec: str
        nreact: str
        speclist: List[str]
        nnz: str = None

    @dataclass
    class InfoContent:
        method: str
        device: str

    @dataclass
    class UserdataContent:
        var: List[str]
        user_var: List[str]

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
        var: List[str] = None

    def __init__(self, netinfo: object, *args, **kwargs) -> None:

        self._env = None
        self._macros = None
        self._info = None
        self._userdata = None
        self._ode = None

    @abstractmethod
    def _prepare_contents(self, netinfo: object, *args, **kwargs) -> None:
        raise NotImplementedError

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

        # self.render_constants(prefix=prefix, save=save)
        self.render_macros(prefix=prefix, save=save)
        self.render_userdata(prefix=prefix, save=save)
        self.render_ode(prefix=prefix, save=save)
        self.render_naunet(prefix=prefix, save=save)

    # This function should not be used outside console commands
    def render_cmake(self, prefix: str = "./", version: str = None) -> None:

        template_names = [
            "CMakeLists.txt.j2",
            "src/CMakeLists.txt.j2",
        ]
        for name in template_names:
            template = self._env.get_template(name)
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
        self, prefix: str = "./", name: str = None, save: bool = True
    ) -> None:
        if not name and save:
            name = "naunet_constants.h"

        template = self._env.get_template("include/naunet_constants.h.j2")
        self._render(
            template, prefix, name, save, constants=self._constants, info=self._info
        )

    def render_macros(
        self, prefix: str = "./", name: str = None, save: bool = True
    ) -> None:
        if not name and save:
            name = "naunet_macros.h"

        template = self._env.get_template("include/naunet_macros.h.j2")
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
            template = self._env.get_template("include/naunet.h.j2")
            self._render(
                template,
                headerprefix,
                headername,
                save,
                info=self._info,
                header=headername,
            )

        template = self._env.get_template("src/naunet.cpp.j2")
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
            template = self._env.get_template("include/naunet_ode.h.j2")
            self._render(
                template, headerprefix, headername, save, ode=self._ode, info=self._info
            )

        template = self._env.get_template("src/naunet_ode.cpp.j2")
        self._render(template, prefix, name, save, ode=self._ode, info=self._info)

    def render_userdata(
        self, prefix: str = "./", name: str = None, save: bool = True
    ) -> None:
        if not name and save:
            name = "naunet_userdata.h"

        template = self._env.get_template("include/naunet_userdata.h.j2")
        self._render(template, prefix, name, save, userdata=self._userdata)


class CVodeTemplateLoader(TemplateLoader):
    def __init__(
        self, netinfo: object, method: str, device: str, *args, **kwargs
    ) -> None:

        super().__init__(netinfo)

        rpath = "templates/cvode"
        loader = PackageLoader("naunet", rpath)
        self._env = Environment(loader=loader)
        self._env.globals.update(zip=zip)
        self._env.trim_blocks = True
        # env.lstrip_blocks = True
        self._env.rstrip_blocks = True

        self._prepare_contents(netinfo, method, device)

    def _prepare_contents(self, netinfo: object, method: str, device: str) -> None:
        n_spec = netinfo.n_spec
        n_react = netinfo.n_react
        species = netinfo.species
        reactions = netinfo.reactions
        databases = netinfo.databases

        self._info = self.InfoContent(method, device)

        nspec = f"#define NSPECIES {n_spec}"
        nreact = f"#define NREACTIONS {n_react}"
        speclist = [f"#define IDX_{x.alias} {i}" for i, x in enumerate(species)]

        self._macros = self.MacrosContent(nspec, nreact, speclist)

        var = [f"double {v};" for db in databases for v in db.variables.values()]
        user_var = []

        self._userdata = self.UserdataContent(var, user_var)

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

            rsym = [y[idx] for idx in rspecidx]
            rsym_mul = "*".join(rsym)
            for specidx in rspecidx:
                rhs[specidx] += f" - {rates[rl]}*{rsym_mul}"
            for specidx in pspecidx:
                rhs[specidx] += f" + {rates[rl]}*{rsym_mul}"

            for specidx in rspecidx:
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
        jac = [
            f"IJth(jmatrix, {idx//n_spec}, {idx%n_spec}) = {j};"
            for idx, j in enumerate(jacrhs)
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

        var = [
            f"realtype {v} = u_data->{v};"
            for db in databases
            for v in db.variables.values()
        ]

        for db in databases:
            if len(db.user_var) > 0:
                var.extend(db.user_var)

        self._ode = self.ODEContent(
            method,
            device,
            rateeqns,
            fex,
            jac,
            spjacrptr=spjacrptr,
            spjaccval=spjaccval,
            spjacdata=spjacdata,
            var=var,
        )


class ODEIntTemplateLoader(TemplateLoader):
    def __init__(
        self, netinfo: object, method: str, device: str, *args, **kwargs
    ) -> None:

        super().__init__(netinfo)

        rpath = "templates/odeint"
        loader = PackageLoader("naunet", rpath)
        self._env = Environment(loader=loader)
        self._env.globals.update(zip=zip)
        self._env.trim_blocks = True
        # env.lstrip_blocks = True
        self._env.rstrip_blocks = True

        self._prepare_contents(netinfo, method, device)

    def _prepare_contents(self, netinfo: object, method: str, device: str) -> None:
        n_spec = netinfo.n_spec
        n_react = netinfo.n_react
        species = netinfo.species
        reactions = netinfo.reactions
        databases = netinfo.databases

        self._info = self.InfoContent(method, device)

        nspec = f"#define NSPECIES {n_spec}"
        nreact = f"#define NREACTIONS {n_react}"
        speclist = [f"#define IDX_{x.alias} {i}" for i, x in enumerate(species)]

        self._macros = self.MacrosContent(nspec, nreact, speclist)

        var = [f"double {v};" for db in databases for v in db.variables.values()]
        user_var = []

        self._userdata = self.UserdataContent(var, user_var)

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

            rsym = [y[idx] for idx in rspecidx]
            rsym_mul = "*".join(rsym)
            for specidx in rspecidx:
                rhs[specidx] += f" - {rates[rl]}*{rsym_mul}"
            for specidx in pspecidx:
                rhs[specidx] += f" + {rates[rl]}*{rsym_mul}"

            for specidx in rspecidx:
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
        jac = [f"j({idx//n_spec}, {idx%n_spec}) = {j};" for idx, j in enumerate(jacrhs)]

        var = [
            f"double {v} = u_data->{v};"
            for db in databases
            for v in db.variables.values()
        ]

        for db in databases:
            if len(db.user_var) > 0:
                var.extend(db.user_var)

        self._ode = self.ODEContent(
            method,
            device,
            rateeqns,
            fex,
            jac,
            var=var,
        )
