from __future__ import annotations
import logging
from pathlib import Path
from typing import Type
from tqdm import tqdm
from .patchmaker import PatchMaker
from .templateloader import NetworkInfo, TemplateLoader
from .species import Species
from .reactions.reaction import Reaction
from .reactions.kidareaction import KIDAReaction
from .reactions.kromereaction import KROMEReaction
from .reactions.leedsreaction import LEEDSReaction
from .reactions.uclchemreaction import UCLCHEMReaction
from naunet.reactions.umistreaction import UMISTReaction
from .dusts.dust import Dust
from .dusts.unidust import UniDust
from .dusts.rr07dust import RR07Dust
from .thermalprocess import ThermalProcess, get_allowed_cooling, get_allowed_heating


logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()


supported_dust_model = {
    "uniform": UniDust,
    "RR07": RR07Dust,
}

supported_reaction_class = {
    "kida": KIDAReaction,
    "umist": UMISTReaction,
    "leeds": LEEDSReaction,
    "krome": KROMEReaction,
    "uclchem": UCLCHEMReaction,
}


def _reaction_factory(react_string: str, format: str) -> Reaction:
    """
    Factory of reactions

    Args:
        react_string (str): the reactions read from file
        format (str): the source of the reaction, use to interpret string format

    Returns:
        Reaction: a reaction object
    """

    initializer = supported_reaction_class.get(format)
    react_string = initializer.preprocessing(react_string)
    if react_string:
        return initializer(react_string)
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
        dusttype: str[str] = None,
    ) -> None:

        self.format_list = set()
        self.reaction_list = []
        self._reactants = set()
        self._products = set()
        self._skipped_reactions = []
        self._info = None
        self._patchmaker = None
        self._templateloader = None

        dust_model = supported_dust_model.get(dusttype)  # dust model class
        # Instantiate a dust model
        self._dust = dust_model() if dust_model else None
        self._allowed_species = allowed_species.copy() if allowed_species else []
        self._required_species = required_species.copy() if required_species else []
        self._allowed_heating = None
        self._allowed_cooling = None
        self._shielding = shielding if shielding else {}
        self._heating_names = heating.copy() if heating else []
        self._cooling_names = cooling.copy() if cooling else []

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
                    for fname, db in zip(filelist, fileformats):
                        self.add_reaction_from_file(fname, db)

            elif isinstance(fileformats, str):

                for fname in filelist:
                    self.add_reaction_from_file(fname, fileformats)

            else:
                raise RuntimeError(f"Unknown format: {fileformats}")

        elif isinstance(filelist, str) or isinstance(filelist, Path):

            filelist = str(filelist)

            if isinstance(fileformats, list):
                if len(fileformats) != 1:
                    raise RuntimeError(
                        "Sizes of input files and sources are mismatching."
                    )
                else:
                    self.add_reaction_from_file(filelist, fileformats[0])

            elif isinstance(fileformats, str):
                self.add_reaction_from_file(filelist, fileformats)

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

        if not isinstance(reaction, Reaction):
            rclass.finalize()

        # reset network information if content is changed
        self._info = None

    def add_reaction_from_file(self, filename: str, format: str) -> None:
        """Add reactions into network from file

        Args:
            filename (str): the file to be read
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
            "The following reactions are duplicate:\n{}".format(
                "\n".join([str(x) for x in dupes])
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

    @property
    def info(self):
        if self._info:
            return self._info

        rqdsp = set([Species(s) for s in self._required_species])
        speclist = sorted(self._reactants | self._products | rqdsp)

        heating = [self.allowed_heating.get(h) for h in self._heating_names]
        cooling = [self.allowed_cooling.get(c) for c in self._cooling_names]

        rclasses = [supported_reaction_class.get(fmt) for fmt in self.format_list]
        self._info = NetworkInfo(
            len(speclist),
            len(self.reaction_list),
            speclist,
            self.reaction_list,
            rclasses,
            heating=heating,
            cooling=cooling,
            dust=self.dust,
            shielding=self._shielding,
        )

        logger.info(
            "{} species in the network: {}".format(
                self._info.n_spec, ", ".join([x.name for x in self._info.species])
            )
        )

        logger.info(
            "Skipped reactions: {}".format(
                "\n".join([str(x) for x in self._skipped_reactions])
            )
        )

        return self._info

    def patchmaker(self, target: str, device: str, *args, **kwargs):

        if self._info and self._patchmaker:
            return self._patchmaker

        self._patchmaker = PatchMaker(self.info, target, device, *args, **kwargs)
        return self._patchmaker

    @property
    def products(self):
        return self._products

    @property
    def reactants(self):
        return self._reactants

    def templateloader(
        self,
        solver: str,
        method: str,
        device: str,
        ratemodifier: dict[int, str] = None,
        odemodifier: list[str] = None,
    ) -> TemplateLoader:

        if self._info and self._templateloader:
            return self._templateloader

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

    def write(self, filename: str, format: str = "") -> None:
        with open(filename, "w") as outf:
            for reac in self.reaction_list:
                outf.write(f"{reac:{format}}\n")
