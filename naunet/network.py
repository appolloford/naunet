from __future__ import annotations
import os
import shutil
import logging
from pathlib import Path
from typing import Type
from tqdm import tqdm
from importlib.metadata import version
from .patchmaker import PatchMaker
from .templateloader import NetworkInfo, TemplateLoader
from .species import Species
from .reactions.reaction import Reaction
from .reactions.kidareaction import KIDAReaction
from .reactions.kromereaction import KROMEReaction
from .reactions.leedsreaction import LEEDSReaction
from .reactions.uclchemreaction import UCLCHEMReaction
from .reactions.umistreaction import UMISTReaction
from .reactions.converter import ExpressionConverter
from .dusts.dust import Dust
from .dusts.hh93dust import HH93Dust
from .dusts.rr07dust import RR07Dust
from .thermalprocess import ThermalProcess, get_allowed_cooling, get_allowed_heating
from .configuration import Configuration


logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()


supported_dust_model = {
    "hh93": HH93Dust,
    "rr07": RR07Dust,
}

supported_reaction_class = {
    "naunet": Reaction,
    "kida": KIDAReaction,
    "umist": UMISTReaction,
    "leeds": LEEDSReaction,
    "krome": KROMEReaction,
    "uclchem": UCLCHEMReaction,
}


def _dust_factory(model: str = "none", **kwargs) -> Dust:

    dustmodel = supported_dust_model.get(model)
    if not dustmodel:
        raise ValueError(f"Unknown dust model: {model}")

    return dustmodel(**kwargs)


def _reaction_factory(react_string: str, format: str) -> Reaction:
    """
    Factory of reactions

    Args:
        react_string (str): the reactions read from file
        format (str): the source of the reaction, use to interpret string format

    Returns:
        Reaction: a reaction object
    """

    initializer = supported_reaction_class.get(format, Reaction)
    react_string = initializer.preprocessing(react_string)
    if react_string:
        return initializer(react_string=react_string)
    return None


def define_reaction(name: str):
    """
    Decorator for users to add customized reaction class

    Args:
        name (str): name of the class / format of the reaction
    """

    def insert_class(reactcls: Type[Reaction]):
        supported_reaction_class.update({name: reactcls})

    return insert_class


def define_dust(name: str):
    """
    Decorator for users to add customized dust model

    Args:
        name (str): name of the class / model of the dust
    """

    def insert_class(dustcls: Type[Dust]):
        supported_dust_model.update({name: dustcls})

    return insert_class


class Network:

    _rateconverter = ExpressionConverter("C")

    def __init__(
        self,
        reactions: list[Reaction] = None,
        filelist: str | list[str] | Path | list[Path] = None,
        fileformats: str | list[str] = None,
        allowed_species: list[str] = None,
        required_species: list[str] = None,
        heating: list[str] = None,
        cooling: list[str] = None,
        shielding: dict[str, str] = None,
        dustmodel: str = "none",
        dustparams: dict = None,
    ) -> None:

        self.format_list = set()
        self.reaction_list = []
        self._reactants = set()
        self._products = set()
        self._skipped_reactions = []
        self._info = None
        self._patchmaker = None
        self._templateloader = None

        dustparams = dustparams if dustparams else {}
        dustparams = {**dustparams, "model": dustmodel}
        self._dust = None if dustmodel == "none" else _dust_factory(**dustparams)
        self._allowed_species = allowed_species.copy() if allowed_species else []
        self._required_species = required_species.copy() if required_species else []
        self._allowed_heating = None
        self._allowed_cooling = None
        self._shielding = shielding if shielding else {}
        self._heating_names = heating.copy() if heating else []
        self._cooling_names = cooling.copy() if cooling else []

        self._reactconsts = {}
        self._reactvaris = {}
        self._reactlocvars = []

        if allowed_species and required_species:

            conflict = [sp for sp in required_species if sp not in allowed_species]

            if conflict:

                raise RuntimeError(
                    """All required species must exist in the allowed species list. 
                    Otherwise leave one of them to be "None"."""
                )

        if reactions:

            if filelist or fileformats:
                logger.warning("Initialize from reactions. Files are skipped!")

            for reac in reactions:
                self.add_reaction(reac)

        # check the lists of files and formats are matching
        # and add them into the network if possible
        elif isinstance(filelist, list):

            filelist = [str(f) for f in filelist]

            if isinstance(fileformats, list):

                if len(filelist) != len(fileformats):
                    raise RuntimeError(
                        "Sizes of input files and sources are mismatching."
                    )
                else:
                    for fname, fmt in zip(filelist, fileformats):
                        if fname and fmt:
                            self.add_reaction_from_file(fname, fmt)

            elif isinstance(fileformats, str):

                for fname in filelist:
                    if fname and fileformats:
                        self.add_reaction_from_file(fname, fileformats)

            else:
                raise RuntimeError(f"Unknown format: {fileformats}")

        elif isinstance(filelist, str) or isinstance(filelist, Path):

            fname = str(filelist)

            if isinstance(fileformats, list):
                if len(fileformats) != 1:
                    raise RuntimeError(
                        "Sizes of input files and sources are mismatching."
                    )
                else:
                    fmt = fileformats[0]
                    if fname and fmt:
                        self.add_reaction_from_file(fname, fmt)

            elif isinstance(fileformats, str):
                if fname and fileformats:
                    self.add_reaction_from_file(fname, fileformats)

            else:
                raise RuntimeError(f"Unknown format: {fileformats}")

        elif filelist is None:
            pass

        else:
            raise TypeError(f"Unknown type of filelist {type(filelist)}")

    def _add_reaction(self, reaction: Reaction | tuple[str, str]) -> list:

        if not isinstance(reaction, Reaction):
            reaction = _reaction_factory(*reaction)

        # return empty set for updating if it is a fake react_string
        if not reaction:
            return set()

        if self._allowed_species:
            if not all(
                [
                    rp.name in self._allowed_species
                    for rp in reaction.reactants + reaction.products
                ]
            ):
                self._skipped_reactions.append(reaction)
                return set()

        self.reaction_list.append(reaction)
        new_reactants = set(reaction.reactants).difference(self._reactants)
        new_products = set(reaction.products).difference(self._products)
        self._reactants.update(new_reactants)
        self._products.update(new_products)
        # if len(self.reaction_list) % 100 == 0:
        #     print("Processing: {} reactions...".format(len(self.reaction_list)))
        return new_reactants | new_products

    def _consts(self) -> dict[str, str]:

        # consts in dust
        dust = self.dust
        dconsts = {f"{c:<15}": f"{cv}" for c, cv in dust.consts.items()} if dust else {}

        # consts in heating/cooling
        thermproc = []
        thermproc.extend([self.allowed_heating.get(h) for h in self._heating_names])
        thermproc.extend([self.allowed_cooling.get(c) for c in self._cooling_names])
        tconsts = {f"{c:<15}": f"{cv}" for p in thermproc for c, cv in p.consts.items()}

        # create consts of binding energies
        rqdsp = set([Species(s) for s in self._required_species])
        speclist = sorted(self._reactants | self._products | rqdsp)
        ebs = {f"eb_{s.alias:<12}": f"{s.eb}" for s in speclist if s.is_surface}

        consts = {**self._reactconsts, **dconsts, **tconsts, **ebs}

        return consts

    def _locvars(self) -> dict[str, str]:

        dust = self.dust
        dlocvars = dust.locvars if dust else []

        thermproc = []
        thermproc.extend([self.allowed_heating.get(h) for h in self._heating_names])
        thermproc.extend([self.allowed_cooling.get(c) for c in self._cooling_names])
        tlocvars = [v for p in thermproc for v in p.locvars]

        locvars = [*self._reactlocvars, *dlocvars, *tlocvars]
        return locvars

    def _varis(self) -> dict[str, str]:

        # varis in dust
        dust = self.dust
        dvaris = {f"{var}": val for var, val in dust.varis.items()} if dust else {}

        # varis in heating/cooling
        thermproc = []
        thermproc.extend([self.allowed_heating.get(h) for h in self._heating_names])
        thermproc.extend([self.allowed_cooling.get(c) for c in self._cooling_names])
        tvaris = {f"{var}": val for p in thermproc for var, val in p.varis.items()}

        varis = {**self._reactvaris, **dvaris, **tvaris}
        return varis

    def add_reaction(self, reaction: Reaction | tuple[str, str]) -> None:
        """Add a reaction into the network

        Args:
            reaction (Reaction | tuple[str, str]): the reaction to be added, either an
                instance of Reaction or a tuple of (reaction_string, format)

        Raises:
            RuntimeError: if the format is unknown when trying to create reaction
                instance from reaction string.
        """

        format = reaction.format if isinstance(reaction, Reaction) else reaction[1]

        self.format_list.update({format} if format else {})
        rclass = supported_reaction_class.get(format)

        if not isinstance(reaction, Reaction):
            # create reaction instance from string
            # change some global settings or class attibutes if needed
            if rclass:
                rclass.initialize()
            else:
                raise RuntimeError(f"Unknown format: {format}")

        new_species = self._add_reaction(reaction)
        logger.info("New species are added: {}".format(new_species))

        self._reactconsts.update(
            {f"{c:<15}": f"{cv}" for c, cv in rclass.consts.items()}
        )
        self._reactvaris.update({f"{var}": val for var, val in rclass.varis.items()})
        self._reactlocvars.extend(
            [v for v in rclass.locvars if v not in self._reactlocvars]
        )

        if not isinstance(reaction, Reaction):
            rclass.finalize()

        # reset network information if content is changed
        self._info = None

    def add_reaction_from_file(self, filename: str | Path, format: str) -> None:
        """Add reactions into network from file

        Args:
            filename (str | Path): the file to be read
            format (str): the format will be used to parse the file

        Raises:
            RuntimeError: if the format is unknown
        """

        self.format_list.update({format})
        new_species = set()

        # change some global settings or class attibutes if needed
        rclass = supported_reaction_class.get(format)
        if rclass:
            rclass.initialize()
        else:
            raise RuntimeError(f"Unknown format: {format}")

        with open(filename, "r") as networkfile:
            for _, line in enumerate(
                tqdm(networkfile.readlines(), desc="Reading File...")
            ):
                new_species.update(self._add_reaction((line, format)))

            # print("New species: \n{}".format("\n".join(str(x) for x in new_species)))

        self._reactconsts.update(
            {f"{c:<15}": f"{cv}" for c, cv in rclass.consts.items()}
        )
        self._reactvaris.update({f"{var}": val for var, val in rclass.varis.items()})
        self._reactlocvars.extend(
            [v for v in rclass.locvars if v not in self._reactlocvars]
        )

        rclass.finalize()

        self._info = None

    @property
    def allowed_cooling(self) -> dict[str, ThermalProcess]:
        """
        Based on the network species to get a dict of allowed cooling processes from
        all support cooling models.

        Returns:
            dict[str, ThermalProcess]: allowed cooling processes.
            Dictionary of {name: <cooling process>}
        """
        if self._allowed_cooling:
            return self._allowed_cooling

        rqdsp = set([Species(s) for s in self._required_species])
        speclist = sorted(self._reactants | self._products | rqdsp)

        self._allowed_cooling = get_allowed_cooling(speclist)
        return self._allowed_cooling

    @property
    def allowed_heating(self) -> dict[str, ThermalProcess]:
        """
        Based on the network species to get a dict of allowed heating processes from
        all support heating models.

        Returns:
            dict[str, ThermalProcess]: allowed heating processes.
                                       Dictionary of {name: <heating process>}
        """
        if self._allowed_heating:
            return self._allowed_heating

        rqdsp = set([Species(s) for s in self._required_species])
        speclist = sorted(self._reactants | self._products | rqdsp)

        self._allowed_heating = get_allowed_heating(speclist)
        return self._allowed_heating

    @property
    def allowed_species(self):
        return self._allowed_species

    @allowed_species.setter
    def allowed_species(self, speclist: list):
        self._allowed_species = speclist.copy()

    @property
    def dust(self):
        return self._dust

    @dust.setter
    def dust(self, dust: Dust) -> None:
        if self._dust and self.reaction_list:
            logger.warning(
                """
                Dust model existed. The existed reactions may 
                not be consistent with newly added reactions
                """
            )
        self._dust = dust

    def check_duplicate_reaction(self, full_check: bool = True):

        seen = {}
        dupes = []

        check_list = (
            self.reaction_list
            if full_check
            else [f"{react:short}" for react in self.reaction_list]
        )

        for x in check_list:
            if x not in seen:
                seen[x] = 1
            else:
                if seen[x] == 1:
                    dupes.append(x)
                seen[x] += 1

        logger.info(
            "The following {} reactions are duplicate:\n{}".format(
                len(dupes), "\n".join([str(x) for x in dupes])
            )
        )
        return dupes

    def check_source_sink(self):
        source = self._reactants.difference(self._products)
        sink = self._products.difference(self._reactants)
        if len(source) == 0 and len(sink) == 0:
            print("Found no source or sink")
        elif len(source) != 0:
            print("Found sources: ", source)
        elif len(sink) != 0:
            print("Found sinks: ", sink)

    def where(
        self,
        reaction: Reaction = None,
        species: Species | str = None,
        mode: str = "all",
    ) -> list[int]:

        indices = []

        if reaction:
            if mode == "all":
                indices = [
                    idx
                    for idx, reac in enumerate(self.reaction_list)
                    if reac == reaction
                ]

            elif mode == "short":
                indices = [
                    idx
                    for idx, reac in enumerate(self.reaction_list)
                    if f"{reac:short}" == f"{reaction:short}"
                ]

            else:
                raise RuntimeError("Unknown mode: {mode}")

        elif species:
            species = species if isinstance(species, Species) else Species(species)

            if mode == "reactant":
                indices = [
                    idx
                    for idx, reac in enumerate(self.reaction_list)
                    if species in reac.reactants
                ]

            elif mode == "product":

                indices = [
                    idx
                    for idx, reac in enumerate(self.reaction_list)
                    if species in reac.products
                ]

            elif mode == "all":

                indices = [
                    idx
                    for idx, reac in enumerate(self.reaction_list)
                    if species in reac.reactants + reac.products
                ]

            else:
                raise RuntimeError("Unknown mode: {mode}")

        return indices

    def export(
        self,
        solver: str = "cvode",
        method: str = "dense",
        device: str = "cpu",
        ratemodifier: dict[int, str] = None,
        odemodifier: list[str] = None,
        prefix: str | Path = "network_export",
        overwrite: bool = False,
    ) -> None:
        tl = self.templateloader(solver, method, device, ratemodifier, odemodifier)

        prefix = Path.cwd() / Path(prefix)
        if not os.path.exists(prefix):
            os.mkdir(prefix)

        elif not overwrite:
            logger.warning("Export directory exists! Stop exporting!")
            return

        reaction_file = prefix / "reactions.naunet"
        if os.path.exists(reaction_file) and not overwrite:
            logger.warning("Reaction file exists! Stop exporting!")
            return

        self.write(reaction_file, "naunet")

        config_file = prefix / "naunet_config.toml"
        if os.path.exists(config_file) and not overwrite:
            logger.warning("Config file exists! Stop exporting!")
            return

        binding = {s.name: s.eb for s in self.info.species if s.is_surface}
        yields = {s.name: s.photon_yield() for s in self.info.species if s.is_surface}

        config = Configuration(
            "network_export",
            description="Exported_network",
            element=Species.known_elements(),
            pseudo_element=Species.known_pseudoelements(),
            allowed_species=self._allowed_species,
            required_species=self._required_species,
            binding_energy=binding,
            photon_yield=yields,
            network=["reactions.naunet"],
            format=["naunet"],
            heating=self._heating_names,
            cooling=self._cooling_names,
            shielding=self._shielding,
            dustmodel=self.dust.model if self.dust else "none",
            dustspecies=Species.dust_species(),
            rate_modifier=ratemodifier,
            ode_modifier=odemodifier,
            solver=solver,
            device=device,
            method=method,
            networkinfo=self.info,
        )

        content = config.content

        with open(config_file, "w", encoding="utf-8") as outf:
            outf.write(content)

        for subdir in ["include", "src", "tests"]:
            subprefix = prefix / subdir

            if os.path.exists(subprefix):

                if not overwrite:
                    logger.warning("Files exist in export directory! Stop exporting!")
                    return

            else:
                os.mkdir(subprefix)

        header_prefix = prefix / "include"
        source_prefix = prefix / "src"
        test_prefix = prefix / "tests"

        tl.render_constants(prefix=source_prefix, headerprefix=header_prefix)
        tl.render_macros(prefix=header_prefix)
        tl.render_naunet(prefix=source_prefix, headerprefix=header_prefix)
        tl.render_ode(prefix=source_prefix, headerprefix=header_prefix)
        tl.render_physics(prefix=source_prefix, headerprefix=header_prefix)
        tl.render_utilities(prefix=source_prefix, headerprefix=header_prefix)
        tl.render_data(prefix=header_prefix)
        tl.render_tests(prefix=test_prefix)
        tl.render_cmake(prefix=prefix, version=version("naunet"))

        pkgpath = Path(__file__).parent
        incfile = pkgpath / "templates/base/cpp/include/naunet_timer.h"
        dest = header_prefix / "naunet_timer.h"
        shutil.copyfile(incfile, dest)

        demo = prefix / "demo.ipynb"
        if not demo.exists():
            shutil.copyfile(pkgpath / "templates/base/demo.ipynb", demo)

    @property
    def info(self):
        if self._info:
            return self._info

        rqdsp = set([Species(s) for s in self._required_species])
        speclist = sorted(self._reactants | self._products | rqdsp)

        heating = [self.allowed_heating.get(h) for h in self._heating_names]
        cooling = [self.allowed_cooling.get(c) for c in self._cooling_names]

        # rclasses = [supported_reaction_class.get(fmt) for fmt in self.format_list]
        self._info = NetworkInfo(
            speclist,
            self.reaction_list,
            heating=heating,
            cooling=cooling,
            dust=self.dust,
            shielding=self._shielding,
            consts=self._consts(),
            varis=self._varis(),
            locvars=self._locvars(),
        )

        nspec = len(self._info.species)
        logger.info(
            "{} species in the network: {}".format(
                nspec, ", ".join([x.name for x in self._info.species])
            )
        )

        logger.info(
            "Skipped reactions: {}".format(
                "\n".join([str(x) for x in self._skipped_reactions])
            )
        )

        return self._info

    def patchmaker(self, target: str, device: str, *args, **kwargs):

        self._patchmaker = PatchMaker(self.info, target, device, *args, **kwargs)
        return self._patchmaker

    @property
    def products(self):
        return self._products

    @property
    def reactants(self):
        return self._reactants

    def reindex(self) -> None:
        for idx, reac in enumerate(self.reaction_list):
            reac.idxfromfile = idx

    def templateloader(
        self,
        solver: str,
        method: str,
        device: str,
        ratemodifier: dict[int, str] = None,
        odemodifier: list[str] = None,
    ) -> TemplateLoader:

        reactindices = [reac.idxfromfile for reac in self.reaction_list]
        if all([idx == -1 for idx in reactindices]):
            self.reindex()
            logger.warning(
                "Reactions have no index information, reindex with the joining order."
            )
        elif any([idx == -1 for idx in reactindices]) and ratemodifier:
            logger.warning(
                "Some reaction has not set index. The rate modifier will not work on them."
            )
        # TODO: repeat indices

        self._templateloader = TemplateLoader(
            self.info, solver, method, device, ratemodifier, odemodifier
        )
        return self._templateloader

    def to_code(
        self,
        solver: str = "cvode",
        method: str = "dense",
        device: str = "cpu",
        ratemodifier: dict[int, str] = None,
        odemodifier: list[str] = None,
        prefix: str = "./",
    ):

        tl = self.templateloader(solver, method, device, ratemodifier, odemodifier)
        tl.render(prefix=prefix, save=True)

    def write(self, filename: str | Path, format: str = "") -> None:

        with open(filename, "w") as outf:

            if format == "krome":
                outf.write("@format:idx,r,r,r,p,p,p,p,p,tmin,tmax,rate\n")

            for reac in self.reaction_list:
                outf.write(f"{reac:{format}}")

                if format == "krome":
                    self._rateconverter.read(reac.rateexpr(self.dust))
                    outf.write(f",{self._rateconverter:fortran}\n")

                else:
                    outf.write(f"\n")
