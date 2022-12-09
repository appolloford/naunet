from __future__ import annotations
import logging
import os
import shutil
from pathlib import Path
from typing import Type
from tqdm import tqdm
from .templateloader import TemplateLoader
from .species import Species
from .reactions import Reaction, builtin_reaction_format
from .reactions.converter import ExpressionConverter
from .reactiontype import ReactionType
from .grains import Grain, builtin_grain_model
from .thermalprocess import ThermalProcess, get_allowed_cooling, get_allowed_heating
from .configuration import Configuration


logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()


supported_grain_model = {cls.model: cls for cls in builtin_grain_model}

supported_reaction_class = {cls.format: cls for cls in builtin_reaction_format}


def _grain_factory(model: str, **kwargs) -> Grain:

    modelcls = supported_grain_model.get(model)
    if model and not modelcls:
        raise ValueError(f"Unknown grain model: {model}")
    else:
        return None if modelcls is None else modelcls(**kwargs)


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

        if not isinstance(name, str):
            raise TypeError("Reaction format name must be string")

        if name in supported_reaction_class.keys():
            raise RuntimeError("Format name has existed, try another name")

        reactcls.format = name
        supported_reaction_class.update({name: reactcls})
        return reactcls

    return insert_class


def define_grain(name: str):
    """
    Decorator for users to add customized grain model
    """

    def insert_class(graincls: Type[Grain]):

        if not isinstance(name, str):
            raise TypeError("Grain model name must be string")

        if name in supported_grain_model.keys():
            raise RuntimeError("Model name has existed, try another name")

        graincls.model = name
        supported_grain_model.update({name: graincls})
        return graincls

    return insert_class


class Network:
    """
    Class of chemical network

    This class generated network instances. The instance can be created from
    text network files with a pre-defined reaction format or from a list of
    reaction instances. It is also allowed to create an empty instance and
    add reactions afterward.

    Examples:
        1. ```network = Network(filelist=<filename>, fileformats="kida")```
        2. ```network = Network([Reaction(...)])```
        3. ```network = Network()
              network.add(Reaction(...))
           ```

    Attributes:
        reaction_list (list[Reaction]): the reactions in the network

    """

    _rateconverter = ExpressionConverter("C")

    def __init__(
        self,
        reactions: list[Reaction] = None,
        filelist: str | list[str] | Path | list[Path] = None,
        fileformats: str | list[str] = None,
        elements: list[str] = None,
        pseudo_elements: list[str] = None,
        allowed_species: list[str] = None,
        allowed_species_kwargs: dict[str, str] = None,
        required_species: list[str] = None,
        required_species_kwargs: dict[str, str] = None,
        heating: list[str] = None,
        cooling: list[str] = None,
        shielding: dict[str, str] = None,
        grain_model: str = "",
        rate_modifier: dict[int, str] = None,
        ode_modifier: dict[str, dict[str, list[str | list[str]]]] = None,
    ) -> None:

        self.reaction_list = []
        self._reactants = set()
        self._products = set()
        self._skipped_reactions = []

        self._elements = elements or []
        self._pseudo_elements = pseudo_elements or []

        if self._elements or self._pseudo_elements:
            Species.set_known_elements(self._elements)
            Species.set_known_pseudoelements(self._pseudo_elements)

        allowed_species = allowed_species or []
        required_species = required_species or []
        allowed_species_kwargs = allowed_species_kwargs or {}
        required_species_kwargs = required_species_kwargs or {}
        self._allowed_species = [
            Species(s, **allowed_species_kwargs) for s in allowed_species
        ]
        self._required_species = [
            Species(s, **required_species_kwargs) for s in required_species
        ]
        self._allowed_species_kwargs = allowed_species_kwargs.copy()
        self._required_species_kwargs = required_species_kwargs.copy()
        self._allowed_heating = None
        self._allowed_cooling = None
        self._shielding = shielding if shielding else {}
        self._heating_names = heating.copy() if heating else []
        self._cooling_names = cooling.copy() if cooling else []
        self._grain_model = grain_model
        self._rate_modifier = rate_modifier.copy() if rate_modifier else {}
        self._ode_modifier = ode_modifier.copy() if ode_modifier else {}

        if self._allowed_species and self._required_species:

            conflict = [
                sp for sp in self._required_species if sp not in self._allowed_species
            ]

            if conflict:

                raise RuntimeError(
                    "All required species must exist in the allowed species list."
                    "Otherwise leave one of them to be 'None'."
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

    def __contains__(self, reac: Reaction) -> bool:
        if not isinstance(reac, Reaction):
            return NotImplemented
        return reac in self.reaction_list

    def _add_reaction(
        self, reaction: Reaction | tuple[str, str]
    ) -> tuple[set[Species], set[Species], Reaction]:

        if not isinstance(reaction, Reaction):
            reaction = _reaction_factory(*reaction)

        # return empty set for updating if it is a fake react_string
        if not reaction:
            return set(), set(), None

        if self._allowed_species:
            if not all(
                [
                    rp in self._allowed_species
                    for rp in reaction.reactants + reaction.products
                ]
            ):
                self._skipped_reactions.append(reaction)
                return set(), set(), None

        self.reaction_list.append(reaction)
        new_reactants = set(reaction.reactants).difference(self._reactants)
        new_products = set(reaction.products).difference(self._products)
        self._reactants.update(new_reactants)
        self._products.update(new_products)
        # if len(self.reaction_list) % 100 == 0:
        #     print("Processing: {} reactions...".format(len(self.reaction_list)))
        return new_reactants, new_products, reaction

    def add_reaction(self, reaction: Reaction | tuple[str, str]) -> None:
        """Add a reaction into the network

        Args:
            reaction (Reaction | tuple[str, str]): the reaction to be added, either an
                instance of Reaction or a tuple of (reaction_string, format)

        Raises:
            RuntimeError: if the format is unknown when trying to create reaction
                instance from reaction string.
        """

        if self._elements or self._pseudo_elements:
            Species.set_known_elements(self._elements)
            Species.set_known_pseudoelements(self._pseudo_elements)

        format = reaction.format if isinstance(reaction, Reaction) else reaction[1]

        rclass = supported_reaction_class.get(format)
        if not isinstance(reaction, Reaction):
            # create reaction instance from string
            # change some global settings or class attibutes if needed
            if rclass:
                rclass.initialize()
            else:
                raise RuntimeError(f"Unknown format: {format}")

        new_reactants, new_products, reactinst = self._add_reaction(reaction)
        if new_reactants:
            logger.info(f"New reactants are added: {new_reactants}")
        if new_products:
            logger.info(f"New products are added: {new_products}")

        if not isinstance(reaction, Reaction):
            rclass.finalize()

    def add_reaction_from_file(self, filename: str | Path, format: str) -> None:
        """Add reactions into network from file

        Args:
            filename (str | Path): the file to be read
            format (str): the format will be used to parse the file

        Raises:
            RuntimeError: if the format is unknown
        """

        if self._elements or self._pseudo_elements:
            Species.set_known_elements(self._elements)
            Species.set_known_pseudoelements(self._pseudo_elements)

        new_reactants = set()
        new_products = set()

        # change some global settings or class attibutes if needed
        rclass = supported_reaction_class.get(format)
        if rclass:
            rclass.initialize()
        else:
            raise RuntimeError(f"Unknown format: {format}")

        with open(filename, "r") as networkfile:
            for ln, line in enumerate(
                tqdm(networkfile.readlines(), desc="Reading File...")
            ):
                try:
                    reac, prod, reactinst = self._add_reaction((line, format))
                    new_reactants.update(reac)
                    new_products.update(prod)
                except Exception as e:
                    logger.error(f"Get error in line {ln}: {line}")
                    raise e

        rclass.finalize()

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

        speclist = sorted(
            self._reactants | self._products | set(self._required_species)
        )

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

        speclist = sorted(
            self._reactants | self._products | set(self._required_species)
        )

        self._allowed_heating = get_allowed_heating(speclist)
        return self._allowed_heating

    @property
    def allowed_species(self) -> list[str]:
        """
        The names of allowed species, which are used to constraint system to
        be solved.

        Returns:
            list[str]: name of allowed species
        """
        return [s.name for s in self._allowed_species]

    @allowed_species.setter
    def allowed_species(self, speclist: list[str]):
        if self._elements or self._pseudo_elements:
            Species.set_known_elements(self._elements)
            Species.set_known_pseudoelements(self._pseudo_elements)

        self._allowed_species = [
            Species(s, **self._allowed_species_kwargs) for s in speclist
        ]

        # examine all reactions again
        recorded_reactions = self.reaction_list + self._skipped_reactions
        self._reactants.clear()
        self._products.clear()
        self.reaction_list = []
        self._skipped_reactions = []

        for reaction in recorded_reactions:
            self.add_reaction(reaction)

    @property
    def cooling(self) -> list[ThermalProcess]:
        """
        The included cooling processes

        Returns:
            list[ThermalProcess]: cooling processes to be used
        """
        return [self.allowed_cooling.get(c) for c in self._cooling_names]

    @property
    def elements(self) -> list[Species]:
        return [spec for spec in self.species if spec.is_atom]

    def export(
        self,
        solver: str = "cvode",
        method: str = "dense",
        device: str = "cpu",
        ratemodifier: dict[int, str] = None,
        odemodifier: dict[str, dict[str, list[str | list[str]]]] = None,
        prefix: str | Path = "network_export",
        overwrite: bool = False,
    ) -> None:
        tl = TemplateLoader(solver, method, device)

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

        binding = {s.name: s.eb for s in self.species if s.is_surface}
        yields = {s.name: s.photon_yield() for s in self.species if s.is_surface}

        config = Configuration(
            "network_export",
            description="Exported_network",
            element=self._elements,
            pseudo_element=self._pseudo_elements,
            allowed_species=self.allowed_species,
            required_species=self.required_species,
            binding_energy=binding,
            photon_yield=yields,
            network=["reactions.naunet"],
            format=["naunet"],
            heating=self._heating_names,
            cooling=self._cooling_names,
            shielding=self._shielding,
            grain_model=self._grain_model,
            rate_modifier=ratemodifier,
            ode_modifier=odemodifier,
            solver=solver,
            device=device,
            method=method,
            instance=self,
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

        tl.render(self, path=prefix)
        tl.render_tests(path=prefix)

        pkgpath = Path(__file__).parent

        demo = prefix / "demo.ipynb"
        if not demo.exists():
            shutil.copyfile(pkgpath / "templates/base/demo.ipynb", demo)

    def find_duplicate_reaction(self, mode: str = "") -> list[tuple[int, Reaction]]:
        """
        Find the duplicate reactions in the network. The default behaviour checks the
        reactants, products, temperature ranges, and reaction type. Only checking the
        reactants and products if mode == "short".

        Args:
            mode (str, optional): the method of checking duplicates. Defaults to "".

        Returns:
            list[tuple[int, Reaction]]: the duplicate reactions and their indices
        """

        seen = {}
        dupes = []

        check_list = (
            [f"{react:{mode}}" for react in self.reaction_list]
            if mode
            else self.reaction_list
        )

        for idx, chk in enumerate(check_list):
            if chk not in seen:
                seen[chk] = 1
            else:
                if seen[chk] >= 1:
                    dupes.append((idx, chk))
                seen[chk] += 1

        return dupes

    def find_source_sink(self) -> tuple[set[Species], set[Species]]:
        """
        Find the source/sink species in the network

        Returns:
            tuple[set[Species], set[Species]]: the tuple of (source, sink), source and
                sink are set of Species
        """
        source = self._reactants.difference(self._products)
        sink = self._products.difference(self._reactants)

        return source, sink

    @property
    def grains(self) -> list[Grain]:
        grain_groups = self.grain_groups

        gspec = [s for s in self.species if s.is_grain]

        if not gspec:
            grains = [
                _grain_factory(
                    self._grain_model,
                    group=group,
                )
                for group in grain_groups
            ]

        else:
            grains = [
                _grain_factory(
                    self._grain_model,
                    species=[g for g in gspec if g.grain_group == group],
                    group=group,
                )
                for group in grain_groups
            ]

        grains = [g for g in grains if g is not None]

        return grains

    @property
    def grain_groups(self) -> list[int]:
        species = self._required_species + list(self._reactants | self._products)
        grain_groups = set([s.grain_group for s in species if s.is_grain])
        surface_groups = set([s.surface_group for s in species if s.is_surface])
        if grain_groups == surface_groups:
            return grain_groups
        elif grain_groups and not surface_groups:
            logging.warning(f"Found grains but no surface reaction is involved")
            return grain_groups
        elif all(g in surface_groups for g in grain_groups):
            logging.warning(
                f"Not enough groups of grains {grain_groups} found to match "
                f"with the surface species groups {surface_groups}."
                f"Surface species groups ({len(surface_groups)}) will be used"
            )
            return surface_groups
        else:
            raise RuntimeError(
                f"Grain groups {grain_groups} cannot match"
                f"with surface groups {surface_groups}"
            )

    @property
    def grain_model(self) -> str:
        return self._grain_model

    @grain_model.setter
    def grain_model(self, model: str) -> None:
        self._grain_model = model

    @property
    def heating(self) -> list[ThermalProcess]:
        """
        The included heating processes

        Returns:
            list[ThermalProcess]: the heating processed to be used
        """
        return [self.allowed_heating.get(h) for h in self._heating_names]

    @property
    def ode_modifier(self) -> dict[str, dict[str, list[str | list[str]]]]:
        return self._ode_modifier

    @ode_modifier.setter
    def ode_modifier(self, omod: dict[str, dict[str, list[str | list[str]]]]) -> None:
        self._ode_modifier = omod.copy()

    @property
    def products(self):
        return self._products

    @property
    def rate_modifier(self) -> dict[int, str]:
        return self._rate_modifier

    @rate_modifier.setter
    def rate_modifier(self, rmod: dict[int, str]) -> None:
        self._rate_modifier = rmod.copy()

    @property
    def reactants(self):
        return self._reactants

    @property
    def reactions(self):
        """
        The reactions in the network. Used for rendering codes. If no reaction
        exists, a dummy reaction is auto-filled. To access the real reactions,
        check the `reaction_list`.

        Returns:
            list[Reaction]: the list of reactions
        """
        return self.reaction_list or [Reaction(reaction_type=ReactionType.DUMMY)]

    def reindex(self) -> None:
        """
        Reindex reactions in the network. The index of reactions will be reset
        to the order in the interior list.
        """
        for idx, reac in enumerate(self.reaction_list):
            reac.idxfromfile = idx

    def remove_reaction(self, reaction: int | Reaction | list[int | Reaction]) -> None:
        """
        Remove one or multiple reaction.

        Args:
            reaction (int | Reaction | list[int] | list[Reaction]): the reactions to be
                removed or their indices of reactions in the network

        Raises:
            TypeError: if the argument is not the type of int | Reaction or
                list[int | Reaction]
        """

        if isinstance(reaction, int):
            self.reaction_list.pop(reaction)

        elif isinstance(reaction, list) and all(isinstance(r, int) for r in reaction):
            self.reaction_list = [
                r for idx, r in enumerate(self.reaction_list) if idx not in reaction
            ]

        elif isinstance(reaction, Reaction):
            self.reaction_list = [r for r in self.reaction_list if r != reaction]

        elif isinstance(reaction, list) and all(
            isinstance(r, Reaction) for r in reaction
        ):
            self.reaction_list = [r for r in self.reaction_list if r not in reaction]

        else:
            raise TypeError

    @property
    def required_species(self) -> list[str]:
        """
        The names of required species, which does not involve in reactions
        but may be required by heating or cooling

        Returns:
            list[str]: required species
        """
        return [s.name for s in self._required_species]

    @required_species.setter
    def required_species(self, speclist: list[str]):
        if self._elements or self._pseudo_elements:
            Species.set_known_elements(self._elements)
            Species.set_known_pseudoelements(self._pseudo_elements)

        self._required_species = [
            Species(s, **self._required_species_kwargs) for s in speclist
        ]

    @property
    def shielding(self) -> dict[str, str]:
        return self._shielding

    @property
    def sources(self) -> list[str]:
        return set([r.source for r in self.reaction_list])

    @property
    def species(self) -> list[Species]:
        """
        Species exists in the network, including the species found from
        reactions and the required species.

        Returns:
            list[Species]: species in the network
        """

        speclist = sorted(
            self._reactants | self._products | set(self._required_species)
        )

        connection = {sp: set() for sp in speclist}
        for reac in self.reaction_list:
            reacrp = reac.reactants + reac.products
            for rp in reacrp:
                connection[rp].update(reacrp)

        speclist = sorted(speclist, key=lambda x: (len(connection[x]), x))

        return speclist

    def to_code(
        self,
        solver: str = "cvode",
        method: str = "dense",
        device: str = "cpu",
        path: str = "./",
    ) -> None:

        tl = TemplateLoader(solver, method, device)
        tl.render(self, path=path, save=True)

    def where_reaction(self, reaction: Reaction, mode: str = "all") -> list[int]:
        """
        Find the index of a reaction

        Args:
            reaction (Reaction): target reaction
            mode (str, optional): the way to compare the equality. Defaults to "all".

        Raises:
            RuntimeError: if mode is unkown

        Returns:
            list[int]: the index(es) of the reaction
        """

        indices = []

        if mode == "all":
            indices = [
                idx for idx, reac in enumerate(self.reaction_list) if reac == reaction
            ]

        elif mode == "short":
            indices = [
                idx
                for idx, reac in enumerate(self.reaction_list)
                if f"{reac:short}" == f"{reaction:short}"
            ]

        else:
            raise RuntimeError("Unknown mode: {mode}")

        return indices

    def where_species(self, species: Species | str, mode: str = "all") -> list[int]:
        """
        Find the index of reactions involving the species. Use `mode` to select find
        in `reactants`, `products`, or `all`.

        Args:
            species (Species | str): target species or its name
            mode (str, optional): where to find species. Defaults to "all".

        Raises:
            RuntimeError: if mode is unkown

        Returns:
            list[int]: the index(es) of the reactions
        """

        indices = []

        if self._elements or self._pseudo_elements:
            Species.set_known_elements(self._elements)
            Species.set_known_pseudoelements(self._pseudo_elements)

        species = (
            species
            if isinstance(species, Species)
            else Species(species, **self._allowed_species_kwargs)
        )

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
                idx for idx, reac in enumerate(self.reaction_list) if species in reac
            ]

        else:
            raise RuntimeError("Unknown mode: {mode}")

        return indices

    def write(self, filename: str | Path, format: str = "") -> None:

        with open(filename, "w") as outf:

            if format == "krome":
                outf.write("@format:idx,r,r,r,p,p,p,p,p,tmin,tmax,rate\n")

            for reac in self.reaction_list:
                outf.write(f"{reac:{format}}")

                if format == "krome":
                    grain_dict = (
                        {g.group: g for g in self.grains} if self.grains else {}
                    )
                    self._rateconverter.read(
                        reac.rateexpr(grain_dict.get(reac.grain_group))
                    )
                    outf.write(f",{self._rateconverter:fortran}\n")

                else:
                    outf.write(f"\n")
