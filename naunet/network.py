from __future__ import annotations
import logging
from dataclasses import dataclass
from typing import Type
from tqdm import tqdm
from .patchmaker import PatchMaker
from .templateloader import TemplateLoader
from .species import Species
from .reactions.reaction import Reaction
from .reactions.kidareaction import KIDAReaction
from .reactions.leedsreaction import LEEDSReaction
from .reactions.kromereaction import KROMEReaction
from .reactions.uclchemreaction import UCLCHEMReaction
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
    "leeds": LEEDSReaction,
    "krome": KROMEReaction,
    "uclchem": UCLCHEMReaction,
}


def _reaction_factory(react_string: str, database: str, dust: Dust = None) -> Reaction:
    """
    Factory of reactions

    Args:
        react_string (str): the reactions read from file
        database (str): the source of the reaction, use to interpret string format
        dust (Dust, optional): the associated dust model, if any. Defaults to None.

    Returns:
        Reaction: a reaction object
    """

    initializer = supported_reaction_class.get(database)
    react_string = initializer.preprocessing(react_string)
    if react_string:
        return initializer(react_string, dust=dust)
    return None


def define_reaction(name: str):
    """
    Decorator for users to add customized reaction class

    Args:
        name (str): name of the class / source database of the reaction
    """

    def insert_class(reactcls: Type[Reaction]):
        supported_reaction_class.update({name: reactcls})

    return insert_class


def define_dust(name: str):
    """
    Decorator for users to add customized dust model

    Args:
        name (str): name of the class / source database of the reaction
    """

    def insert_class(dustcls: Type[Dust]):
        supported_dust_model.update({name: dustcls})

    return insert_class


class Network:
    @dataclass
    class Info:
        n_spec: int
        n_react: int
        species: list[Species]
        reactions: list[Reaction]
        databases: list[Type[Reaction]]
        heating: list[str] = None
        cooling: list[str] = None
        dust: Dust = None
        shielding: dict = None
        odemodifier: list[str] = None
        ratemodifier: list[str] = None

    def __init__(
        self,
        filelist: str | list[str] = None,
        filesources: str | list[str] = None,
        allowed_species: list[str] = None,
        required_species: list[str] = None,
        heating: list[str] = None,
        cooling: list[str] = None,
        shielding: dict[str, str] = None,
        dusttype: str[str] = None,
        ratemodifier: dict[int, str] = None,
    ) -> None:

        self.database_list = set()
        self.reaction_list = []
        self.reactants_in_network = set()
        self.products_in_network = set()
        self.ode_modifier = []
        self._rate_modifier = ratemodifier.copy() if ratemodifier else {}
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

        # check the lists of files and databases are matching
        # and add them into the network if possible
        if isinstance(filelist, list):

            if isinstance(filesources, list):

                if len(filelist) != len(filesources):
                    raise RuntimeError(
                        "Sizes of input files and sources are mismatching."
                    )
                else:
                    for fname, db in zip(filelist, filesources):
                        self.add_reaction_from_file(fname, db)

            elif isinstance(filesources, str):

                for fname in filelist:
                    self.add_reaction_from_file(fname, filesources)

            else:
                raise RuntimeError(f"Unknown format: {filesources}")

        elif isinstance(filelist, str):

            if isinstance(filesources, list):
                if len(filesources) != 1:
                    raise RuntimeError(
                        "Sizes of input files and sources are mismatching."
                    )
                else:
                    self.add_reaction_from_file(filelist, filesources[0])

            elif isinstance(filesources, str):
                self.add_reaction_from_file(filelist, filesources)

            else:
                raise RuntimeError(f"Unknown format: {filesources}")

    def _add_reaction(self, react_string: str, database: str) -> list:

        reaction = _reaction_factory(react_string, database, self._dust)

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
        new_reactants = set(reaction.reactants).difference(self.reactants_in_network)
        new_products = set(reaction.products).difference(self.products_in_network)
        self.reactants_in_network.update(new_reactants)
        self.products_in_network.update(new_products)
        # if len(self.reaction_list) % 100 == 0:
        #     print("Processing: {} reactions...".format(len(self.reaction_list)))
        return new_reactants | new_products

    def add_reaction(self, react_string: str, database: str) -> None:

        self.database_list.update({database})

        # change some global settings or class attibutes if needed
        rclass = supported_reaction_class.get(database)
        if rclass:
            rclass.initialize()
        else:
            raise RuntimeError(f"Unknown format: {database}")

        new_species = self._add_reaction(react_string, database)
        logger.info("New species are added: {}".format(new_species))

        rclass.finalize()

        # reset network information if content is changed
        self._info = None

    def add_reaction_from_file(self, filename: str, database: str) -> None:

        if not filename:
            logger.critical("No file assigned!")

        if not database:
            logger.critical(
                """
                Try to read in file but database is not assigned. 
                Try again by "add_reaction_from_file"
                """
            )

        self.database_list.update({database})
        new_species = set()

        # change some global settings or class attibutes if needed
        rclass = supported_reaction_class.get(database)
        if rclass:
            rclass.initialize()
        else:
            raise RuntimeError(f"Unknown format: {database}")

        with open(filename, "r") as networkfile:
            for _, line in enumerate(
                tqdm(networkfile.readlines(), desc="Reading File...")
            ):
                new_species.update(self._add_reaction(line, database))

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
        speclist = sorted(self.reactants_in_network | self.products_in_network | rqdsp)

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
        speclist = sorted(self.reactants_in_network | self.products_in_network | rqdsp)

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

    @property
    def info(self):
        if self._info:
            return self._info

        rqdsp = set([Species(s) for s in self._required_species])
        speclist = sorted(self.reactants_in_network | self.products_in_network | rqdsp)

        heating = [self.allowed_heating.get(h) for h in self._heating_names]
        cooling = [self.allowed_cooling.get(c) for c in self._cooling_names]

        databaselist = [supported_reaction_class.get(db) for db in self.database_list]
        self._info = self.Info(
            len(speclist),
            len(self.reaction_list),
            speclist,
            self.reaction_list,
            databaselist,
            heating=heating,
            cooling=cooling,
            dust=self.dust,
            shielding=self._shielding,
            odemodifier=self.ode_modifier,
            ratemodifier=self._rate_modifier,
        )

        logger.info(
            "{} species in the network: {}".format(
                self._info.n_spec, ", ".join([x.name for x in self._info.species])
            )
        )

        logger.info(
            "Skipped reactions: {}".format(
                "\n".join([repr(x) for x in self._skipped_reactions])
            )
        )

        return self._info

    def patchmaker(self, target: str, device: str, *args, **kwargs):

        if self._info and self._patchmaker:
            return self._patchmaker

        self._patchmaker = PatchMaker(self.info, target, device, *args, **kwargs)
        return self._patchmaker

    @property
    def ratemodifier(self) -> dict[int, str]:
        return self._rate_modifier

    @ratemodifier.setter
    def ratemodifier(self, ratemodifier: dict[int, str]) -> None:
        self._rate_modifier = ratemodifier.copy()

    def templateloader(self, solver: str, method: str, device: str):

        if self._info and self._templateloader:
            return self._templateloader

        self._templateloader = TemplateLoader(self.info, solver, method, device)
        return self._templateloader

    def to_code(
        self,
        solver: str = "cvode",
        method: str = "dense",
        device: str = "cpu",
        prefix: str = "./",
    ):

        tl = self.templateloader(solver, method, device)
        tl.render(prefix=prefix, save=True)
