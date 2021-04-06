import logging
from tqdm import tqdm
from . import settings
from .templateloader import ODEIntTemplateLoader, TemplateLoader, CVodeTemplateLoader
from .reactions.reaction import Reaction
from .reactions.kidareaction import KIDAReaction
from .reactions.leedsreaction import LEEDSReaction
from .reactions.kromereaction import KROMEReaction


logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()


supported_reaction_class = {
    "kida": KIDAReaction,
    "leeds": LEEDSReaction,
    "krome": KROMEReaction,
}

supported_template_loader = {
    "cvode": CVodeTemplateLoader,
    "odeint": ODEIntTemplateLoader,
}


def reaction_factory(react_string: str, database: str) -> Reaction:
    initializer = supported_reaction_class.get(database)
    react_string = initializer.preprocessing(react_string)
    if react_string:
        return initializer(react_string, format)
    return None


class Network:
    class Info:
        def __init__(self, species: list, reactions: list, databases: list) -> None:
            self.n_spec = len(species)
            self.n_react = len(reactions)
            self.species = species
            self.reactions = reactions
            self.databases = databases

    def __init__(
        self,
        fname: str = None,
        database: str = None,
        species: list = None,
    ) -> None:

        self.database_list = set()
        self.reaction_list = []
        self.reactants_in_network = set()
        self.products_in_network = set()
        self._allowed_species = None
        self._skipped_reactions = []
        self._info = None
        self._templateloader = None

        if species:
            self._allowed_species = species

        if fname and database:
            self.add_reaction_from_file(fname, database)

    def _add_reaction(self, react_string: str, database: str) -> list:
        reaction = reaction_factory(react_string, database)

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
        new_species = self._add_reaction(react_string, database)
        logger.info("New species are added: {}".format(new_species))

        # reset network information if content is changed
        self._info = None

    def add_reaction_from_file(self, filename: str, database: str) -> None:
        if not filename:
            logger.critical("No file assigned!")
        if not database:
            logger.critical(
                'Try to read in file but database is not assigned. Try again by "add_reaction_from_file"'
            )
        self.database_list.update({database})
        new_species = set()
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
        self._info = self.Info(speclist, self.reaction_list, databaselist)
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

    def templateloader(self, solver: str, method: str, device: str):

        if self._info and self._templateloader:
            return self._templateloader

        tl = supported_template_loader.get(solver)
        self._templateloader = tl(self.info, method, device)
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
