from __future__ import annotations
import logging
import re
from . import chemistrydata


class Species:
    """
    Class of chemical species

    This class generated species instances from the name of species. The
    name of species are like H2O, C2H6, etc. It parses the species name and
    find the substring matches with the recorded element/pseudo-element
    names. The list of element names can be set by the users (see methods
    below). Otherwise the default lists of element/pseudo-element are used.

    Notes:
        1. Pseudo-element here is used to say the special symbols exist in
        species but not elements in periodic table, isotopes, charges. E.g.
        o(ortho-), p(para-), CR (cosmic-ray)
        2. For species start with `grain_symbol` or `surface_prefix`. They
        will be identified as grain and surface species (species stick
        on grains).
            - Species("GC", surface_prefix="G") is equivalent to Species("#C")
        3. Charges must be +/- symbol append in the end. E.g. He++, Si+++

    Attributes:
        default_elements (list): The default element names used to parse
            species names.
        default_pseudoelement (list): The default pseudo-element names
            used to parse species names.
        element_count (dict[str, int]): The count of each element in the
            known species list.

    """

    default_elements = [
        "e",
        "E",
        "H",
        "D",
        "He",
        "C",
        "N",
        "O",
        "F",
        "Na",
        "Mg",
        "Al",
        "Si",
        "P",
        "S",
        "Cl",
        "Ar",
        "Ca",
        "Fe",
        "Ni",
    ]
    default_pseudoelements = [
        "CR",
        "CRP",
        "XRAY",
        "Photon",
        "PHOTON",
        "CRPHOT",
        "X",
        "M",
        "p",
        "o",
        "m",
        "c-",
        "l-",
        r"\*",
        "g",
    ]

    _known_elements = []
    _known_pseudoelements = []
    _replacement = {}

    def __init__(
        self,
        name: str,
        grain_symbol: str = "GRAIN",
        surface_prefix: str = "#",
        bulk_prefix: str = "@",
    ) -> None:
        """
        Initializes Species with species name

        Args:
            name (str): name of the species
        """

        if not isinstance(name, str):
            raise TypeError("Species need a name string to instantiate")

        self.name = name
        self.element_count = dict()
        self._alias = None
        self._binding_energy = None
        self._enthalpy = None
        self._mass = 0.0
        self._massnumber = 0
        self._photon_yield = None
        self._is_grain = False
        self._grain_symbol = grain_symbol
        self._grain_group = None
        self._is_surface = False
        self._surface_prefix = surface_prefix
        self._surface_group = None
        self._bulk_prefix = bulk_prefix

        # Initialize known elements if not set when the fist species is instanciated
        if not Species._known_elements and not Species._known_pseudoelements:
            logging.warning("No assigned known element list. Use default elements")

            Species._known_elements.extend(Species.default_elements)
            Species._known_pseudoelements.extend(Species.default_pseudoelements)

        self._parse_molecule_name(
            Species._known_elements + Species._known_pseudoelements,
            [self._grain_symbol, self._surface_prefix],
        )

    def __copy__(self) -> Species:
        return type(self)(
            self.name,
            self._grain_symbol,
            self._surface_prefix,
            self._bulk_prefix,
        )

    def __eq__(self, o: Species) -> bool:
        if isinstance(o, Species):
            return (
                (self.is_electron and o.is_electron)
                or (
                    self.is_grain
                    and o.is_grain
                    and self.grain_group == o.grain_group
                    and self.charge == o.charge
                )
                or (
                    self.is_surface
                    and o.is_surface
                    and self.surface_group == o.surface_group
                    and self.charge == o.charge
                    and self.basename == o.basename
                )
                or self.name == o.name
            )
            # return (self.is_electron and o.is_electron) or self.name == o.name
        return NotImplemented

    def __hash__(self) -> int:
        return (
            hash("Electron")
            if self.is_electron
            else hash(
                f"{self.basename}"
                f"{self.charge}"
                f"{self.is_grain}"
                f"{self.grain_group}"
                f"{self.is_surface}"
                f"{self.surface_group}"
            )
        )

    def __lt__(self, o) -> bool:
        if isinstance(o, Species):
            return self.name < o.name
        return NotImplemented

    def __repr__(self) -> str:
        return (
            f"Species('{self.name}', "
            f"grain_symbol='{self._grain_symbol}', "
            f"surface_prefix='{self._surface_prefix}', "
            f"bulk_prefix='{self._bulk_prefix}')"
        )

    def __format__(self, spec) -> str:
        return f"{self.name:{spec}}"

    def _add_element_count(self, element: str, count: int):
        """
        Add the count of elements in the species. Used in `_parse_molecule_name`

        Args:
            element (str): The name of the elements
            count (int): the number to be added into count
        """
        if element in self._known_pseudoelements:
            return

        elif element == self._surface_prefix:
            if self._surface_group is not None:
                raise RuntimeError(f"Repeatedly found surface symbol.")
            else:
                self._is_surface = True
                self._surface_group = count
                return

        elif element == self._grain_symbol:
            if self._grain_group is not None:
                raise RuntimeError(
                    f"Repeatedly found grain symbol {self._grain_symbol} in {self.name}."
                )
            else:
                self._is_grain = True
                self._grain_group = count
                count = 1

        elif count <= 0:
            logging.warning(
                f"The number added into {element} count <= 0 in {self.name}. Reset to 1"
            )
            count = 1
        if element in self.element_count.keys():
            self.element_count[element] += count
        else:
            self.element_count[element] = count

    @classmethod
    def _check_elements(cls) -> None:
        """
        Check whether repeated elements exist in `_known_elements` and
        `_known_pseudoelements` lists.
        """
        all_element = cls._known_elements + cls._known_pseudoelements
        dup_element = [ele for ele in all_element if all_element.count(ele) > 1]
        if dup_element:
            logging.error(
                f"Repeated definitions of {dup_element} in elements and pseudoelemts"
            )

    def _parse_molecule_name(self, elements: list[str], symbols: list[str]) -> None:
        """
        Parse the name of the species to get the composition of elements.

        Raises:
            RuntimeError: [description]
        """

        parsename = self.name
        # remove charge symbols
        parsename = re.sub(r"\+*$", "", parsename)
        parsename = re.sub(r"-*$", "", parsename)
        charge = self.name.replace(parsename, "")

        components = sorted(elements + symbols, key=len, reverse=True)

        firstparse = parsename
        matches = []
        for c in components:
            for it in re.finditer(c, firstparse):
                matches.append(it)
                start, end = it.start(), it.end()
                # remove the found items to avoid repeatance (e.g. S in Si)
                firstparse = firstparse[:start] + " " * (end - start) + firstparse[end:]
        matches = sorted(matches, key=lambda x: x.start())

        # to check the number of elements, we need the end of last match and
        # the start of the next one
        starts = [m.start() for m in matches]
        ends = [m.end() for m in matches]
        # make the ends start from 0 and starts end at string length so that
        # we can check the number of the last element. (e.g CO2)
        starts.append(len(parsename))
        ends.insert(0, 0)
        # insert dummy name to match the length of starts and ends
        matchnames = [m.group() for m in matches]
        matchnames.insert(0, "")

        if starts[0] != 0:
            raise RuntimeError(f"{self.name} starts with something unrecognizable")

        for s, e, n in zip(starts, ends, matchnames):
            # if there is replacement, save the element name with the new value
            n = self._replacement.get(n, n)
            if e != s:
                substring = parsename[e:s]
                if substring.isdigit():
                    self._add_element_count(n, int(parsename[e:s]))
                else:
                    raise RuntimeError(
                        f'Unrecongnized name: "{substring}" in "{self.name}"'
                    )
            else:
                if n in symbols:
                    self._add_element_count(n, 0)
                elif n:
                    self._add_element_count(n, 1)

        # replace the original species name if the element name is in the replace target
        if self._replacement:
            newname = ""
            for s, e, n in zip(starts, ends, matchnames):
                nr = self._replacement.get(n, n)
                newname = f"{newname}{nr}{parsename[e:s]}"
            self.name = f"{newname}{charge}"

    @classmethod
    def add_known_elements(cls, elements: list[str]) -> list:
        """
        Add names of elements to the list of known elements

        Args:
            elements (list[str]): names of elements

        Raises:
            TypeError: If argument is not a list

        Returns:
            list: the current list of known elements
        """
        if not isinstance(elements, list):
            raise TypeError(f"{elements} is not a list")

        for ele in elements:
            if ele in cls._known_elements:
                logging.warning(f"{ele} exists in element list, skip!")
            elif ele in cls._known_pseudoelements:
                logging.warning(
                    f"{ele} exists in pseudo element list, move to element list!"
                )
                cls._known_pseudoelements.remove(ele)
                cls._known_elements.append(ele)
            else:
                cls._known_elements.append(ele)
        cls._check_elements()
        return cls._known_elements

    @classmethod
    def add_known_pseudoelements(cls, pelements: list[str]) -> list:
        """
        Add names of elements to the list of known pseudo elements

        Args:
            pelements (list[str]): names of pseudo elements

        Raises:
            TypeError: If argument is not a list

        Returns:
            list: the current list of known pseudo elements
        """
        if not isinstance(pelements, list):
            raise TypeError(f"{pelements} is not a list")

        for ele in pelements:
            if ele in cls._known_pseudoelements:
                logging.warning("{} exists in pseudo element list, skip!".format(ele))
            elif ele in cls._known_elements:
                logging.warning(
                    "{} exists in element list, move to pseudo element list!".format(
                        ele
                    )
                )
                cls._known_elements.remove(ele)
                cls._known_pseudoelements.append(ele)
            else:
                cls._known_pseudoelements.append(ele)
        cls._check_elements()
        return cls._known_pseudoelements

    @property
    def alias(self) -> str:
        """
        Alias of the species.

        Replace surface symbol by `G` and charges by `M`(-), `I`(+). Capital
        element symbols are replaced by standard element symbols. Used in
        indexing species in codes. Can be customized by setter if needed.
        """

        if not self._alias:
            basename = self.basename
            # TODO: The replacement does not guarantee the correctness
            # e.g. CO could be replaced by Co if Co exists in the known element list
            replacement = {
                e.Symbol.upper(): e.Symbol
                for e in chemistrydata.periodic_table + chemistrydata.isotopes_table
                if e.Symbol.upper() in self._known_elements
            }
            for key, value in replacement.items():
                basename = basename.replace(key, value)
            self._alias = "{}{}{}".format(
                "G" if self.is_surface else "",
                basename,
                "I" * (self.charge + 1) if self.charge >= 0 else "M" * abs(self.charge),
            )
        return self._alias

    @alias.setter
    def alias(self, name: str):
        self._alias = name

    @property
    def basename(self) -> str:
        """The species name without surface symbols and charges."""

        basename = self.name
        # remove surface symbol
        if self.is_surface:
            prefix = f"{self._surface_prefix}{self._surface_group or ''}"
            basename = basename.replace(prefix, "")

        # remove charge symbols
        if self.charge != 0:
            basename = re.sub(r"\+*$", "", basename)
            basename = re.sub(r"-*$", "", basename)

        return basename

    @property
    def binding_energy(self) -> float:
        """
        Binding energy (K) of surface species.

        Return the binding energy if this species is at ice-phase. If the
        binding energy was not set before, it will be searched in the UMIST
        2012 binding energy data. Raise error when no value is found.

        Raises:
            RuntimeError: binding energy cannot be found
        """
        if not self.is_surface:
            raise ValueError(f"{self.name} is not ice species, has no binding energy")

        eb = (
            self._binding_energy
            or chemistrydata.user_binding_energy.get(self.name)
            or chemistrydata.rate12_binding_energy.get(self.gasname)
        )

        if not eb:
            raise RuntimeError(f"Cannot find the binding energy of {self.name}")

        return eb

    @binding_energy.setter
    def binding_energy(self, eb: float) -> None:
        self._binding_energy = eb

    # alias function
    eb = binding_energy

    @property
    def charge(self) -> int:
        """
        Total charge of the species. Calculate the number of "+" and "-" in the end
        """
        if self.is_electron:
            return -1

        pcharge = "".join(re.findall(r"\+*$", self.name)).count("+")
        ncharge = "".join(re.findall(r"-*$", self.name)).count("-")
        return pcharge - ncharge

    @property
    def enthalpy(self) -> float:
        """
        The enthalpy (kJ/mol) of formation at 0K (surface species only, from CCCBDB)
        """
        if not self.is_surface:
            return 0.0

        avails = {e.Species: float(e.Enthalpy) for e in chemistrydata.enthalpy_table}

        return self._enthalpy or avails.get(self.gasname, 0.0)

    @enthalpy.setter
    def enthalpy(self, val: float) -> None:
        self._enthalpy = val

    @property
    def gasname(self) -> str:
        """
        Return the name of its gas-phase species if this species is at ice-phase.
        Else return the current name.
        """
        prefix = f"{self._surface_prefix}{self._surface_group or ''}"
        return self.name.replace(prefix, "") if self.is_surface else self.name

    @property
    def grain_group(self) -> int:
        return self._grain_group

    @property
    def is_electron(self) -> bool:
        """
        A conveninent function to check whether the species is electron.

        Returns:
            bool: True if the species is electron
        """
        return self.name.upper() in ["E", "E-"]

    @property
    def is_atom(self) -> bool:
        """
        Check whether the species is the atomic form of an element, including
        grain species (e.g. GRAIN0) if it is provided

        Returns:
            bool: True if the species is an atom
        """
        names, counts = self.element_count.keys(), self.element_count.values()
        return (
            len(names) == 1
            and sum(counts) == 1
            and self.charge == 0
            and not self.is_electron
            and not self.is_surface
        )

    @property
    def is_grain(self) -> bool:
        return self._is_grain

    @property
    def is_surface(self) -> bool:
        """
        Check whether the species is sticking on surface (a surface species)
        """
        return self._is_surface

    # TODO: python 3.9 support classmethod property
    @classmethod
    def known_elements(cls) -> list:
        """
        Returns the current list of known elements

        Returns:
            list: the current list of known elements
        """
        # return the default element list if none of them initialized
        if not cls._known_elements and not cls._known_pseudoelements:
            return cls.default_elements
        else:
            return cls._known_elements

    @classmethod
    def known_pseudoelements(cls) -> list:
        """
        Returns the current list of known pseudo elements

        Returns:
            list: the current list of known pseudo elements
        """
        # return the default pseudo element list if none of them initialized
        if not cls._known_elements and not cls._known_pseudoelements:
            return cls.default_pseudoelements
        else:
            return cls._known_pseudoelements

    @property
    def mass(self) -> float:
        """
        The mass (amu) of the species, estimated by summing the mass of elements
        """
        if self._mass:
            return self._mass

        self._mass = 0.0
        for e in chemistrydata.periodic_table + chemistrydata.isotopes_table:
            self._mass += self.element_count.get(e.Symbol, 0) * float(e.AtomicMass)

        # ? electron mass
        # self._mass -= 0.00054858 * self.charge

        # TODO: PAH, GRAIN, electron
        if self._mass <= 0.0:
            logging.warning(f"{self.name} has zero or negative mass (amu)!")

        return self._mass

    @property
    def massnumber(self) -> float:
        """
        The mass number (neutron + proton) of the species
        """
        if self._massnumber:
            return self._massnumber

        self._massnumber = 0.0
        for e in chemistrydata.periodic_table + chemistrydata.isotopes_table:
            self._massnumber += self.element_count.get(e.Symbol, 0) * (
                float(e.NumberofNeutrons) + float(e.NumberofProtons)
            )

        # TODO: PAH, GRAIN, electron
        # self._massnumber += self.element_count.get("GRAIN", 0) * 1200.0

        if self._massnumber <= 0.0:
            logging.warning(f"{self.name} has zero or negative mass number!")

        return self._massnumber

    # alias function of mass number
    A = massnumber

    @property
    def n_atoms(self) -> int:
        """
        The number of atoms in the species. 0 for electron and 1 for grain.
        """
        if self.is_electron:
            return 0
        names, counts = self.element_count.keys(), self.element_count.values()
        return sum(counts)

    @property
    def photon_yield(self) -> float:
        """
        The photodesorption yield of the species (only for ice-phase species).
        Return 0.0 if no value is found.

        Returns:
            float: photodesorption yield
        """
        if not self.is_surface:
            logging.fatal(f"{self.name} is not ice species! No photodesorption yield")

        return self._photon_yield or chemistrydata.user_photon_yield.get(self.name, 0.0)

    @photon_yield.setter
    def photon_yield(self, phyld: float) -> None:
        self._photon_yield = phyld

    @classmethod
    def remove_known_elements(cls, elements: list) -> list:
        """
        Remove names of elements to the list of known elements

        Args:
            elements (list): names of elements

        Raises:
            TypeError: If argument is not a list

        Returns:
            list: the current list of known elements
        """
        if not isinstance(elements, list):
            raise TypeError(f"{elements} is not a list")

        for ele in elements:
            cls._known_elements.remove(ele)
        cls._check_elements()
        return cls._known_elements

    @classmethod
    def remove_known_pseudoelements(cls, pelements: list) -> list:
        """
        Add names of elements to the list of known pseudo elements

        Args:
            pelements (list): names of pseudo elements

        Raises:
            TypeError: If argument is not a list

        Returns:
            list: the current list of known pseudo elements
        """
        if not isinstance(pelements, list):
            raise TypeError(f"{pelements} is not a list")

        for ele in pelements:
            cls._known_pseudoelements.remove(ele)
        cls._check_elements()
        return cls._known_pseudoelements

    @classmethod
    def reset(cls) -> None:
        """Reset class attributes in Species

        Reset known_elements and known_pseudoelements
        """
        cls._known_elements = []
        cls._known_pseudoelements = []
        cls._replacement = {}

    @classmethod
    def set_known_elements(cls, elements: list) -> None:
        """
        Set the list of known elements to new list

        Args:
            elements (list): names of elements

        Raises:
            TypeError: If argument is not a list

        Returns:
            list: the current list of known elements
        """
        if not isinstance(elements, list):
            raise TypeError(f"{elements} is not a list")

        cls._known_elements.clear()
        cls._known_elements.extend(elements)
        cls._check_elements()

    @classmethod
    def set_known_pseudoelements(cls, pelements: list) -> None:
        """
        Set the list of known pseudo elements to new list

        Args:
            pelements (list): names of pseudo elements

        Raises:
            TypeError: If argument is not a list

        Returns:
            list: the current list of known pseudo elements
        """
        if not isinstance(pelements, list):
            raise TypeError(f"{pelements} is not a list")

        cls._known_pseudoelements.clear()
        cls._known_pseudoelements.extend(pelements)
        cls._check_elements()

    @property
    def surface_group(self) -> int:
        return self._surface_group


def top_abundant_species(
    species_list: list[Species],
    abundances: list[float],
    element: str = None,
    rank: int = -1,
) -> tuple[Species, float]:
    """
    Find the most abundant species, if element is provides, find the main
    reservoirs of the elements.

    Args:
        species_list (list[Species]): list of species, must have the same
            order as abundances
        abundances (list[float]): list of the abundances of species
        element (str, optional): name of target element. Defaults to None.
        rank (int, optional): the top n species to be shown. Defaults to -1
            returns whole list.

    Raises:
        RuntimeError: The element doesn't exist in any species.

    Returns:
        tuple[Species, float]: the most abundant species and their abundances
    """
    sorted_abund = None
    if element == None:
        sorted_abund = sorted(
            zip(species_list, abundances), key=lambda x: x[1], reverse=True
        )[:rank]

    else:
        filtered_abund = list(
            filter(
                lambda x: element in x[0].element_count.keys(),
                zip(species_list, abundances),
            )
        )
        sorted_abund = sorted(
            filtered_abund,
            key=lambda x: x[1] * x[0].element_count[element],
            reverse=True,
        )[:rank]

    if not sorted_abund:
        raise RuntimeError(
            "Undefined results. please check the element exists in the network."
        )

    return sorted_abund
