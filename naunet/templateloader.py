from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from importlib.metadata import version
from pathlib import Path
from tqdm import tqdm
from typing import TYPE_CHECKING
from jinja2 import Template, Environment, PackageLoader

from .species import Species
from .reactions.reaction import Reaction
from .reactiontype import ReactionType
from .thermalprocess import ThermalProcess
from .grains.grain import Grain
from .utilities import _collect_variable_items, _prefix, _suffix, _stmwrap

if TYPE_CHECKING:
    from .network import Network

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
    grains: list[Grain]
    shielding: dict


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
    class Jacobian:
        """
        Analytical Jacobian matrix of the differtial equations. The sparse form
        is save in CSR style.

        Attributes:
            nnz (int): number of non-zero terms
            rows (list[int]): number of elements in a row
            cols (list[int]): column index of elements
            vals (list[str]): non-zero terms
            rhs: all elements, including zero terms
        """

        nrow: int
        nnz: int
        rows: list[int]
        cols: list[int]
        vals: list[str]
        rhs: list[str]

    @dataclass
    class ODEContent:
        """
        The information required by ode expressions
        """

        rateeqns: list[str]
        hrateeqns: list[str]
        crateeqns: list[str]
        fex: list[str]
        jac: TemplateLoader.Jacobian

    @dataclass
    class RenormContent:
        factor: list[str]
        matrix: list[str]

    def __init__(self, solver: str, method: str, device: str) -> None:
        loader = PackageLoader("naunet")
        # self._env = RelativeEnvironment(loader=loader)
        self._env = Environment(loader=loader)
        self._env.globals.update(zip=zip)
        self._env.filters["collect_variable_items"] = _collect_variable_items
        self._env.filters["prefix"] = _prefix
        self._env.filters["suffix"] = _suffix
        self._env.filters["stmwrap"] = _stmwrap
        self._env.trim_blocks = True
        self._env.rstrip_blocks = True

        self._solver = solver
        self._general = self.GeneralInfo(method, device, version("naunet"))

    def _assign_rates(
        self,
        rate_sym: str,
        reactions: list[Reaction | ThermalProcess],
        grains: list[Grain] = None,
    ) -> list[str]:
        # check the temperature range exists
        ltranges = [f"Tgas>={r.temp_min}" if r.temp_min > 0 else "" for r in reactions]
        utranges = [f"Tgas<{r.temp_max}" if r.temp_max > 0 else "" for r in reactions]
        tranges = [
            "".join([lt, " && " if lt and ut else "", ut])
            for lt, ut in zip(ltranges, utranges)
        ]

        grain_dict = {g.group: g for g in grains} if grains else {}

        if grains:
            # chemical reations
            rateexprs = [
                reac.rateexpr(grain_dict.get(reac.grain_group)) for reac in reactions
            ]
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

    def _prepare_ode_content(
        self,
        netinfo: NetworkInfo,
        species_kwargs: dict[str, str] = None,
        rate_modifier: dict[int, str] = None,
        ode_modifier: dict[str, dict[str, list[str | list[str]]]] = None,
    ) -> ODEContent:
        species = netinfo.species
        reactions = netinfo.reactions

        heating = netinfo.heating
        cooling = netinfo.cooling
        grains = netinfo.grains

        has_thermal = True if netinfo.heating or netinfo.cooling else False
        n_spec = len(species)
        n_eqns = max(n_spec + has_thermal, 1)

        species_kwargs = species_kwargs or {}

        rate_sym = "k"
        rateeqns = self._assign_rates(rate_sym, reactions, grains)

        for idx, reac in enumerate(reactions):
            for key, value in rate_modifier.items():
                if key == reac.idxfromfile:
                    logging.warning(f"Overwirte the rate of: `{reac}` with {value}")
                    rateeqns[idx] = f"{rate_sym}[{idx}] = {value};"

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

        # add the modifying term to fex and jac
        # the performance could be bad if fex mismatch with jac
        for sname, expr in ode_modifier.items():
            spec = Species(sname, **species_kwargs)
            sidx = species.index(spec)
            for fact, dep in zip(expr["factors"], expr["reactants"]):
                depspec = [Species(d, **species_kwargs) for d in dep]
                depsym = [f"y[IDX_{d.alias}]" for d in depspec]
                depsym_mul = "*".join(depsym)

                rhs[sidx] += f" + ({fact}) * {depsym_mul}"

                for dspec in depspec:
                    didx = species.index(dspec)
                    depsymcopy = depsym.copy()
                    depsymcopy.remove(y[didx])
                    depsymcopy_mul = "*".join(depsymcopy)

                    term = f" + {'*'.join([f'({fact})', *depsymcopy_mul])}"
                    jacrhs[sidx * n_eqns + didx] += term

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
            rhs[n_spec] = f"(gamma - 1.0) * ( {rhs[n_spec]} ) / kerg / npar"
            for si in range(n_spec):
                jacrhs[n_spec * n_eqns + si] = (
                    "0.0"
                    if jacrhs[n_spec * n_eqns + si] == "0.0"
                    else f"(gamma - 1.0) * ( {jacrhs[n_spec * n_eqns + si]} ) / kerg / npar"
                )

        fex = [f"{l} = {r};" for l, r in zip(lhs, rhs)]

        spjacrptr = []
        spjaccval = []
        spjacdata = []

        nnz = 0

        for row in range(n_eqns):
            spjacrptr.append(nnz)
            for col in range(n_eqns):
                elem = jacrhs[row * n_eqns + col]
                if elem != "0.0":
                    spjaccval.append(col)
                    spjacdata.append(f"{elem}")
                    nnz += 1
        spjacrptr.append(nnz)

        jac = self.Jacobian(n_eqns, nnz, spjacrptr, spjaccval, spjacdata, jacrhs)

        return self.ODEContent(rateeqns, hrateeqns, crateeqns, fex, jac)

    def _prepare_renorm_content(self, netinfo: NetworkInfo) -> RenormContent:
        # get the exact element string
        elements = netinfo.elements
        species = netinfo.species
        elemnames = [next(iter(elem.element_count)) for elem in elements]

        matrix = []
        for iele, einame in enumerate(elemnames):
            for jele, ejname in enumerate(elemnames):
                terms = ["0.0"]
                for ispec, spec in enumerate(species):
                    ci = spec.element_count.get(einame, 0)
                    cj = spec.element_count.get(ejname, 0)
                    if not spec.is_electron and ci and cj:
                        terms.append(
                            f"{(ci * cj * elements[jele].A)} * ab[IDX_{spec.alias}] / {spec.A} / Hnuclei"
                        )
                matrix.append(" + ".join(terms))

        renorm = []
        for spec in species:
            counts = [spec.element_count.get(ename, 0) for ename in elemnames]
            factor = [
                f"{c * elem.A} * rptr[IDX_ELEM_{ename}] / {spec.A}"
                for c, ename, elem in zip(counts, elemnames, elements)
                if c
            ]
            renorm.append(1.0 if spec.is_electron else " + ".join(factor))

        return self.RenormContent(renorm, matrix)

    def _render(
        self,
        template: Template,
        info: NetworkInfo = None,
        ode: ODEContent = None,
        renorm: RenormContent = None,
        save: bool = True,
        path: Path | str = None,
    ) -> None:
        result = template.render(
            general=self._general,
            network=info,
            ode=ode,
            renorm=renorm,
        )
        name = template.name.replace(".j2", "")
        name = name.replace(f"{self._solver}/", "")

        cuda_support = ["constants", "fex", "jac", "physics", "rates", "renorm"]
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
        network: Network,
        templates: list[str] = None,
        save: bool = True,
        path: Path | str = None,
        jac_pattern: bool = False,
    ) -> None:
        templates = templates or self.templates
        solver = self._solver

        reactindices = [reac.idxfromfile for reac in network.reactions]
        if all([idx == -1 for idx in reactindices]):
            network.reindex()
            logging.warning(
                "Reactions have no index information, reindex with the joining order."
            )
        elif any([idx == -1 for idx in reactindices]) and network.rate_modifier:
            logging.warning(
                "Some reaction has not set index. The rate modifier will not work on them."
            )

        info = NetworkInfo(
            network.elements,
            network.species,
            network.reactions or [Reaction(reaction_type=ReactionType.DUMMY)],
            network.heating,
            network.cooling,
            network.grains,
            network.shielding,
        )

        speckws = network._species_kwargs
        rate_modifier = network.rate_modifier
        ode_modifier = network.ode_modifier
        ode = self._prepare_ode_content(info, speckws, rate_modifier, ode_modifier)
        renorm = self._prepare_renorm_content(info)

        for tmplname in templates:
            tmpl = self._env.get_template(f"{solver}/{tmplname}")
            self._render(tmpl, info, ode, renorm, save, path)

        if jac_pattern:
            jacrhs = ode.jac.rhs
            n_eqns = ode.jac.nrow

            pattern = [0 if j == "0.0" else 1 for j in jacrhs]

            rowpattern = []
            for row in range(n_eqns):
                rowdata = pattern[row * n_eqns : (row + 1) * n_eqns]
                rowpattern.append(" ".join(str(e) for e in rowdata))

            pattern = "\n".join(rowpattern)

            with open(path / "jac_pattern.dat", "w") as outf:
                outf.write(pattern)

    @property
    def templates(self) -> list[str]:
        """
        List of all available templates for current solver

        Returns:
            list[str]: list of templates
        """
        solver = self._solver
        return [
            tmpl.replace(f"{solver}/", "")
            for tmpl in self._env.list_templates()
            if tmpl.startswith(solver)
        ]

    def render_tests(self, path: Path | str) -> None:
        path = Path(path)

        testpkgpath = f"templates/base/cpp/tests"
        testenv = Environment(loader=PackageLoader("naunet", testpkgpath))
        tmplnamelist = testenv.list_templates()

        for tmplname in tmplnamelist:
            tmpl = testenv.get_template(tmplname)
            self._render(tmpl, save=True, path=path / "tests")
