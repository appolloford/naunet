from __future__ import annotations

import os
from dataclasses import dataclass
from importlib.metadata import version
from pathlib import Path
from tqdm import tqdm
from jinja2 import Template, Environment, PackageLoader

from .species import Species
from .reactions.reaction import Reaction
from .reactions.reactiontype import ReactionType
from .thermalprocess import ThermalProcess
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
    elements: list[Species]
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
    class GeneralInfo:
        """
        General information
        """

        method: str
        device: str
        version: str

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
        nnz: int
        spjacrptr: list[str] = None
        spjaccval: list[str] = None
        spjacrptrarr: str = None
        spjaccvalarr: str = None
        spjacdata: list[str] = None
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
        self._env.trim_blocks = True
        self._env.rstrip_blocks = True

        self._network_info = netinfo
        self._solver = solver

        self._ode_modifier = odemodifier.copy() if odemodifier else []
        self._rate_modifier = ratemodifier.copy() if ratemodifier else {}

        self._general = self.GeneralInfo(method, device, version("naunet"))
        self._ode = None
        self._physics = None

        self._prepare_contents(netinfo)

    def _assign_rates(
        self,
        rate_sym: str,
        reactions: list[Reaction | ThermalProcess],
        dust: Dust = None,
    ) -> list[str]:

        # check the temperature range exists
        ltranges = [f"Tgas>={r.temp_min}" if r.temp_min > 0 else "" for r in reactions]
        utranges = [f"Tgas<{r.temp_max}" if r.temp_max > 0 else "" for r in reactions]
        tranges = [
            "".join([lt, " && " if lt and ut else "", ut])
            for lt, ut in zip(ltranges, utranges)
        ]

        if dust:
            # chemical reations
            rateexprs = [reac.rateexpr(dust) for reac in reactions]
        else:
            # thermal process
            rateexprs = [reac.rateexpr() for reac in reactions]

        rateassign = [
            "\n".join(
                [
                    f"if ({trange}) {{",
                    f"{rate_sym}[{ridx}] = {rateexpr};",
                    f"}}",
                ]
            )
            if trange
            else f"{rate_sym}[{ridx}] = {rateexpr};"
            for ridx, (trange, rateexpr) in enumerate(zip(tranges, rateexprs))
        ]

        return rateassign

    def _prepare_contents(self, netinfo: NetworkInfo) -> None:

        elements = netinfo.elements
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
        n_eqns = max(n_spec + has_thermal, 1)

        # get the exact element string
        elenames = [next(iter(ele.element_count)) for ele in elements]

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
        )

        rate_sym = "k"
        rateeqns = self._assign_rates(rate_sym, reactions, dust)

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
                rhs[specidx] += f" - {rate_sym}[{rl}]*{rsym_mul}"
            for specidx in pspecidx:
                rhs[specidx] += f" + {rate_sym}[{rl}]*{rsym_mul}"

            # Jacobian
            for specidx in rspecidx:
                # df/dx, remove the dependency for current reactant
                for ri in rspecidx:
                    rsymcopy = rsym.copy()
                    rsymcopy.remove(y[ri])
                    term = f" - {'*'.join([f'{rate_sym}[{rl}]', *rsymcopy])}"
                    jacrhs[specidx * n_eqns + ri] += term
            for specidx in pspecidx:
                for ri in rspecidx:
                    rsymcopy = rsym.copy()
                    rsymcopy.remove(y[ri])
                    term = f" + {'*'.join([f'{rate_sym}[{rl}]', *rsymcopy])}"
                    jacrhs[specidx * n_eqns + ri] += term

            # for specidx, _ in enumerate(species):

            #     # make the diagonal terms nonzero to make sure the matrix reversible
            #     if rhs[specidx] == "0.0":
            #         rhs[specidx] = "naunet_rate_tiny"

        # prepare heating rate expressions
        hrate_sym = "kh"
        hrateeqns = self._assign_rates(hrate_sym, heating)

        # fex/jac of thermal process
        for hidx, h in enumerate(heating):

            rspecidx = [species.index(r) for r in h.reactants]
            rsym = [y[idx] for idx in rspecidx]
            rsym_mul = "*".join(rsym)

            # temperature index is n_spec
            rhs[n_spec] += f" + {hrate_sym}[{hidx}] * {rsym_mul}"

            for ri in rspecidx:
                rsymcopy = rsym.copy()
                rsymcopy.remove(y[ri])
                term = f" + {'*'.join([f'{hrate_sym}[{hidx}]', *rsymcopy])}"
                # only fill the last row of jacobian
                jacrhs[n_spec * n_eqns + ri] += term

        # prepare cooling rate expressions
        crate_sym = "kc"
        crateeqns = self._assign_rates(crate_sym, cooling)

        for cidx, c in enumerate(cooling):

            rspecidx = [species.index(r) for r in c.reactants]
            rsym = [y[idx] for idx in rspecidx]
            rsym_mul = "*".join(rsym)

            # temperature index is n_spec
            rhs[n_spec] += f" - {crate_sym}[{cidx}] * {rsym_mul}"

            for ri in rspecidx:
                rsymcopy = rsym.copy()
                rsymcopy.remove(y[ri])
                term = f" - {'*'.join([f'{crate_sym}[{cidx}]', *rsymcopy])}"
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

        spjacrptrarr = ", ".join(spjacrptrarr)
        spjaccvalarr = ", ".join(spjaccvalarr)

        renormmat = []
        for iele, einame in enumerate(elenames):
            for jele, ejname in enumerate(elenames):
                term = f"IDX_ELEM_{einame}, IDX_ELEM_{ejname}"
                if self._solver == "cvode":
                    term = f"IJth(A, {term}) = 0.0"
                elif self._solver == "odeint":
                    term = f"A({term}) = 0.0"

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
            nnz=nnz,
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
        self,
        template: Template,
        save: bool = True,
        path: Path | str = None,
    ) -> None:

        result = template.render(
            network=self._network_info,
            general=self._general,
            physics=self._physics,
            ode=self._ode,
        )
        name = template.name.replace(".j2", "")
        name = name.replace(f"{self._solver}/", "")

        cuda_support = ["constants", "fex", "jac", "physics", "rates"]
        for substr in cuda_support:
            if substr in name and self._general.device == "gpu":
                name = name.replace("cpp", "cu")

        if save:
            path = Path(path)
            headerpath = path / "include"
            sourcepath = path / "src"

            for p in [path, headerpath, sourcepath]:
                if not p.exists():
                    p.mkdir(parents=True)

            print(path / name)
            with open(path / name, "w") as outf:
                outf.write(result)
        else:
            print(result)

    def render(
        self,
        templates: list[str] = None,
        save: bool = True,
        path: Path | str = None,
    ) -> None:

        templates = templates or self.templates
        solver = self._solver

        for tmplname in templates:
            tmpl = self._env.get_template(f"{solver}/{tmplname}")
            self._render(tmpl, save, path)

    @property
    def templates(self):
        solver = self._solver
        return [
            tmpl.replace(f"{solver}/", "")
            for tmpl in self._env.list_templates()
            if tmpl.startswith(solver)
        ]

    def render_jac_pattern(self, prefix: str = "./") -> None:

        with open(f"{prefix}/jac_pattern.dat", "w") as outf:
            outf.write(self._ode.jacpattern)

    def render_tests(self, path: Path | str) -> None:

        path = Path(path)

        testpkgpath = f"templates/base/cpp/tests"
        testenv = Environment(loader=PackageLoader("naunet", testpkgpath))
        tmplnamelist = testenv.list_templates()

        for tmplname in tmplnamelist:
            tmpl = testenv.get_template(tmplname)
            self._render(tmpl, True, path / "tests")
