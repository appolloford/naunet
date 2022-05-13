from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from tqdm import tqdm
from jinja2 import Template, Environment, PackageLoader

from .species import Species
from .reactions.reaction import Reaction
from .reactions.reactiontype import ReactionType
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
    species: list[Species]
    reactions: list[Reaction]
    heating: list[str]
    cooling: list[str]
    dust: Dust
    shielding: dict
    consts: dict[str, str]
    varis: dict[str, str]
    locvars: list[str]


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

        nreact: int
        nheating: int
        ncooling: int
        elemcidx: list[str]
        speccidx: list[str]
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
        renormmat: list[str] = None
        renorm: list[str] = None
        odemodifier: list[str] = None
        ratemodifier: list[str] = None

    @dataclass
    class PhysicsContent:
        """
        The information required by general physics function
        """

        mantles: str
        mu: str
        gamma: str
        elemabund: list[str]
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
        varis: dict[str, float]
        locvars: list[str]
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

    def _prepare_contents(self, netinfo: NetworkInfo, method: str, device: str) -> None:

        species = netinfo.species
        reactions = netinfo.reactions
        if not reactions:
            reactions.append(Reaction(reaction_type=ReactionType.DUMMY))
            netinfo.varis.update({**Reaction.varis})

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
        n_spec = len(species)
        n_react = len(reactions)
        n_eqns = max(n_spec + has_thermal, 1)

        self._info = self.InfoContent(method, device)

        nheating = len(heating) if heating else 0
        ncooling = len(cooling) if cooling else 0
        # speclist = [f"#define IDX_{x.alias} {i}" for i, x in enumerate(species)]
        speccidx = [f"IDX_{x.alias}" for x in species]
        elements = [spec for spec in species if spec.is_atom]
        # get the exact element string
        elenames = [next(iter(ele.element_count)) for ele in elements]
        elemcidx = [f"IDX_ELEM_{ename}" for ename in elenames]

        consts = netinfo.consts
        varis = netinfo.varis
        locvars = netinfo.locvars

        self._variables = self.VariablesContent(consts, varis, locvars)

        self._macros = self.MacrosContent(
            n_react,
            nheating,
            ncooling,
            elemcidx,
            speccidx,
        )

        mantles = " + ".join(f"y[IDX_{g.alias}]" for g in species if g.is_surface)
        mantles = mantles if mantles else "0.0"

        # TODO: exclude electron, grain?
        density = " + ".join(f"y[IDX_{s.alias}]*{s.massnumber}" for s in species)
        npartile = " + ".join(f"y[IDX_{s.alias}]" for s in species)
        mu = f"({density}) / ({npartile})" if density else "0.0"

        # TODO: different ways to get adiabatic index
        gamma = "5.0 / 3.0"

        elemabund = []
        for iele, einame in enumerate(elenames):
            term = "0.0"
            for ispec, spec in enumerate(species):
                ci = spec.element_count.get(einame, 0)
                if ci:
                    term = f"{term} + {ci:.1f} * y[IDX_{spec.alias}]"
            elemabund.append(term)

        self._physics = self.PhysicsContent(
            mantles,
            mu,
            gamma,
            elemabund,
            h2shielding=netinfo.shielding.get("H2", ""),
            coshielding=netinfo.shielding.get("CO", ""),
            n2shielding=netinfo.shielding.get("N2", ""),
        )

        # prepare reaction rate expressions
        rates = [f"k[{r}]" for r in range(n_react)]
        ltrange = [f"Tgas>={r.temp_min}" if r.temp_min > 0 else "" for r in reactions]
        utrange = [f"Tgas<{r.temp_max}" if r.temp_max > 0 else "" for r in reactions]
        crits = [
            "".join([lt, " && " if lt and ut else "", ut])
            for lt, ut in zip(ltrange, utrange)
        ]
        rateeqns = [
            "\n".join(
                [
                    f"if ({crit}) {{",
                    f"{rate} = {reac.rateexpr(dust)};",
                    f"}}",
                ]
            )
            if crit
            else f"{rate} = {reac.rateexpr(dust)};"
            for crit, rate, reac in zip(crits, rates, reactions)
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
        ltrange = [f"Tgas>={h.temp_min}" if h.temp_min > 0 else "" for h in heating]
        utrange = [f"Tgas<{h.temp_max}" if h.temp_max > 0 else "" for h in heating]
        crits = [
            "".join([lt, " && " if lt and ut else "", ut])
            for lt, ut in zip(ltrange, utrange)
        ]
        hrateeqns = [
            "\n".join(
                [
                    f"if ({crit}) {{",
                    f"{hrate} = {h.rateexpr()};",
                    f"}}",
                ]
            )
            if crit
            else f"{hrate} = {h.rateexpr()};"
            for crit, hrate, h in zip(crits, hrates, heating)
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
        ltrange = [f"Tgas>={c.temp_min}" if c.temp_min > 0 else "" for c in cooling]
        utrange = [f"Tgas<{c.temp_max}" if c.temp_max > 0 else "" for c in cooling]
        crits = [
            "".join([lt, " && " if lt and ut else "", ut])
            for lt, ut in zip(ltrange, utrange)
        ]
        crateeqns = [
            "\n".join(
                [
                    f"if ({crit}) {{",
                    f"{crate} = {c.rateexpr()};",
                    f"}}",
                ]
            )
            if crit
            else f"{crate} = {c.rateexpr()};"
            for crit, crate, c in zip(crits, crates, cooling)
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

        renormmat = []
        for iele, einame in enumerate(elenames):
            for jele, ejname in enumerate(elenames):
                term = f"IJth(A, IDX_ELEM_{einame}, IDX_ELEM_{ejname}) = 0.0"
                for ispec, spec in enumerate(species):
                    ci = spec.element_count.get(einame, 0)
                    cj = spec.element_count.get(ejname, 0)
                    if not spec.iselectron and ci and cj:
                        term = f"{term} + {(ci * cj * elements[jele].A)} * ab[IDX_{spec.alias}] / {spec.A} / Hnuclei"
                renormmat.append(term)

        renorm = []
        for ispec, spec in enumerate(species):
            if spec.iselectron:
                factor = "1.0"
            else:
                factor = "0.0"
                for einame, espec in zip(elenames, elements):
                    ci = spec.element_count.get(einame, 0)
                    if ci:
                        factor = f"{factor} + {(ci * espec.A)} * rptr[IDX_ELEM_{einame}] / {spec.A}"
            renorm.append(f"ab[IDX_{spec.alias}] = ({factor}) * ab[IDX_{spec.alias}]")

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
            renormmat=renormmat,
            renorm=renorm,
            odemodifier=odemodifier,
            ratemodifier=ratemodifier,
        )

    def _render(
        self, template: Template, prefix: str, name: str, save: bool, *args, **kwargs
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
            tname = "base/cpp/include/naunet_constants.h.j2"
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

        tname = "base/cpp/src/naunet_constants.cpp.j2"
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

        tname = "base/cpp/include/naunet_macros.h.j2"
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
            macros=self._macros,
            physics=self._physics,
            variables=self._variables,
            ode=self._ode,
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
            tname = "base/cpp/include/naunet_physics.h.j2"
            template = self._env.get_template(tname)
            self._render(
                template,
                headerprefix,
                headername,
                save,
                physics=self._physics,
                info=self._info,
            )

        tname = "base/cpp/src/naunet_physics.cpp.j2"
        template = self._env.get_template(tname)
        self._render(
            template,
            prefix,
            name,
            save,
            physics=self._physics,
            info=self._info,
            macros=self._macros,
        )

    def render_data(
        self, prefix: str = "./", name: str = None, save: bool = True
    ) -> None:
        if not name and save:
            name = "naunet_data.h"

        tname = "base/cpp/include/naunet_data.h.j2"
        template = self._env.get_template(tname)
        self._render(template, prefix, name, save, variables=self._variables)

    def render_jac_pattern(self, prefix: str = "./") -> None:

        with open(f"{prefix}/jac_pattern.dat", "w") as outf:
            outf.write(self._ode.jacpattern)

    def render_tests(self, prefix) -> None:

        testpkgpath = f"templates/base/cpp/tests"
        testenv = Environment(loader=PackageLoader("naunet", testpkgpath))
        tmplnamelist = testenv.list_templates()

        for tmplname in tmplnamelist:
            tmpl = self._env.get_template(f"base/cpp/tests/{tmplname}")
            name = tmplname.replace(".j2", "")
            self._render(tmpl, prefix, name, True)

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
            tname = "base/cpp/include/naunet_utilities.h.j2"
            template = self._env.get_template(tname)
            self._render(template, headerprefix, headername, save)

        tname = "base/cpp/src/naunet_utilities.cpp.j2"
        template = self._env.get_template(tname)
        self._render(template, prefix, name, save)
