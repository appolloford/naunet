import logging
from tqdm import tqdm
from .patchmaker import PatchMaker
from .templateloader import TemplateLoader
from .reactions.reaction import Reaction
from .reactions.kidareaction import KIDAReaction
from .reactions.leedsreaction import LEEDSReaction
from .reactions.kromereaction import KROMEReaction
from .dusts.dust import Dust
from .dusts.unidust import UniDust


logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()


supported_dust_model = {
    "uniform": UniDust,
}

supported_reaction_class = {
    "kida": KIDAReaction,
    "leeds": LEEDSReaction,
    "krome": KROMEReaction,
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

    def insert_class(reactcls: object):
        supported_reaction_class.update({name: reactcls})

    return insert_class


class Network:
    class Info:
        def __init__(
            self,
            species: list,
            reactions: list,
            databases: list,
            dust: Dust = None,
            shielding: dict = None,
            odemodifier: list = None,
            ratemodifier: list = None,
        ) -> None:

            self.n_spec = len(species)
            self.n_react = len(reactions)
            self.species = species
            self.reactions = reactions
            self.databases = databases
            self.dust = dust
            self.shielding = shielding
            self.odemodifier = odemodifier
            self.ratemodifier = ratemodifier

    # TODO: heating/cooling have not been implemented
    def __init__(
        self,
        species: list = None,
        heating: list = None,
        cooling: list = None,
        shielding: dict = None,
        dusttype: str = None,
    ) -> None:

        self.database_list = set()
        self.reaction_list = []
        self.reactants_in_network = set()
        self.products_in_network = set()
        self.ode_modifier = []
        self.rate_modifier = []
        self._allowed_species = species
        self._skipped_reactions = []
        self._info = None
        self._patchmaker = None
        self._templateloader = None

        # Instantiate a dust model
        self._dust = supported_dust_model.get(dusttype)()
        self._shielding = shielding if shielding else {}

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

        new_species = self._add_reaction(react_string, database)
        logger.info("New species are added: {}".format(new_species))

        # reset network information if content is changed
        self._info = None

    def add_reaction_from_file(self, filename: str, database: str) -> None:

        self.database_list.update({database})
        new_species = set()

        # change some global settings or class attibutes if needed
        rclass = supported_reaction_class.get(database)
        if rclass:
            rclass.initialize()

        with open(filename, "r") as networkfile:
            for _, line in enumerate(
                tqdm(networkfile.readlines(), desc="Reading File...")
            ):
                new_species.update(self._add_reaction(line, database))

            # print("New species: \n{}".format("\n".join(str(x) for x in new_species)))

        self._info = None

    @property
    def allowed_species(self):
        return self._allowed_species

    @allowed_species.setter
    def allowed_species(self, speclist: list):
        self._allowed_species = speclist

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

    def finalize(self):
        for db in list(self.database_list):
            if supported_reaction_class.get(db):
                supported_reaction_class.get(db).finalize()

    @property
    def info(self):
        if self._info:
            return self._info

        speclist = sorted(self.reactants_in_network | self.products_in_network)
        databaselist = [supported_reaction_class.get(db) for db in self.database_list]
        self._info = self.Info(
            speclist,
            self.reaction_list,
            databaselist,
            dust=self.dust,
            shielding=self._shielding,
            odemodifier=self.ode_modifier,
            ratemodifier=self.rate_modifier,
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
