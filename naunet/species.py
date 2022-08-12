from __future__ import annotations
import logging
import re
import csv
import os
from collections import namedtuple
from pathlib import Path
from . import chemistry


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
        2. For surface species (species stick on grains), another special
        symbol are defined in `surface_prefix`. Species starts with the
        special symbol are denoted as surface species.
        3. Charges must be +/- symbol append in the end. E.g. He++, Si+++

    Attributes:
        default_elements (list): The default element names used to parse
            species names.
        default_pseudoelement (list): The default pseudo-element names
            used to parse species names.
        surface_symbol (str): The symbol used to define surface species.
            name (str): Name of the species
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
        "GRAIN",
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
        "*",
        "g",
    ]
    surface_prefix = "#"

    _known_elements = []
    _known_pseudoelements = []
    _periodic_table = []
    _isotopes_table = []
    _enthalpy_table = []
    _dust_species = []

    def __init__(self, name: str) -> None:
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
        self._enthalpy = 0.0
        self._mass = 0.0
        self._massnumber = 0
        self._photon_yield = None
        self._surface_prefix = self.surface_prefix

        # Initialize known elements if not set when the fist species is instanciated
        if not Species._known_elements and not Species._known_pseudoelements:
            logging.warning("No assigned known element list. Use default elements")

            Species._known_elements.extend(Species.default_elements)
            Species._known_pseudoelements.extend(Species.default_pseudoelements)

        chemistry_data_path = Path(chemistry.__file__).parent
        # Create periodic table if not created
        if not Species._periodic_table:
            # Load the periodic table as a list of namedtuple
            # Source: https://gist.github.com/GoodmanSciences/c2dd862cd38f21b0ad36b8f96b4bf1ee
            with open(chemistry_data_path / "periodictable.csv", newline="") as infile:
                # remove comment lines in case
                reader = csv.reader(filter(lambda row: row[0] != "#", infile))
                Element = namedtuple("Element", next(reader))
                Species._periodic_table = list(map(Element._make, reader))

        if not Species._isotopes_table:
            # reference: http://moltensalt.org/references/static/downloads/pdf/stable-isotopes.pdf
            with open(chemistry_data_path / "isotopestable.csv", newline="") as infile:
                # remove comment lines in case
                reader = csv.reader(filter(lambda row: row[0] != "#", infile))
                Isotope = namedtuple("Isotope", next(reader))
                Species._isotopes_table = list(map(Isotope._make, reader))

        if not Species._enthalpy_table:
            # reference: https://cccbdb.nist.gov/hf0k.asp
            with open(chemistry_data_path / "enthalpytable.csv", newline="") as infile:
                # remove comment lines in case
                reader = csv.reader(filter(lambda row: row[0] != "#", infile))
                Enthalpy = namedtuple("Enthalpy", next(reader))
                Species._enthalpy_table = list(map(Enthalpy._make, reader))

        self._parse_molecule_name()

        # create a element count dict in all caps to search in periodic table
        self._allcaps_element_count = {
            key.upper(): value for key, value in self.element_count.items()
        }

    def __copy__(self) -> Species:
        return type(self)(self.name)

    def __eq__(self, o: object) -> bool:
        if isinstance(o, Species):
            return (self.iselectron and o.iselectron) or self.name == o.name
        return NotImplemented

    def __hash__(self) -> int:
        return hash(self.name)

    def __lt__(self, o) -> bool:
        if isinstance(o, Species):
            return self.name < o.name
        return NotImplemented

    def __repr__(self) -> str:
        return f"Species('{self.name}')"

    def __format__(self, spec) -> str:
        return f"{self.name:{spec}}"

    def _add_element_count(self, element: str, count: int):
        """
        Add the count of elements in the species. Used in `_parse_molecule_name`

        Args:
            element (str): The name of the elements
            count (int): the number to be added into count
        """
        if element in Species._known_pseudoelements:
            return
        if count < 0:
            logging.warning(
                "The number added into {} count <= 0 in {}. Reset to 0".format(
                    element, self.name
                )
            )
            count = 0
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

    def _parse_molecule_name(self):
        """
        Parse the name of the species to get the composition of elements.

        Raises:
            RuntimeError: [description]
        """
        element_sorted = sorted(
            Species._known_elements + Species._known_pseudoelements,
            key=len,
            reverse=True,
        )
        element_len = list(map(lambda x: len(x), element_sorted))

        # get the name without surface or charge symbol
        specname = self.basename

        lastelement = None
        while len(specname) > 0:
            # compare the head of the name with the element names
            # remove it after count it
            for ele, l in zip(element_sorted, element_len):
                if specname[0:l] == ele:
                    self._add_element_count(ele, 1)
                    specname = specname[l:]
                    lastelement = ele
                    break
            else:
                # get the leading value and add the lost counts to the last found element
                # without considering the cases like H_2^13CO
                num = re.findall(r"^\d+", specname)
                if len(num) > 0:
                    num = int(num[0])
                    self._add_element_count(lastelement, num - 1)
                    specname = re.sub(r"^\d+", "", specname)
                else:
                    raise RuntimeError(
                        'Unrecongnized name: "{}" in "{}"'.format(specname, self.name)
                    )

    @classmethod
    def add_dust_species(cls, dustspecies: list[str]) -> list[str]:
        if not isinstance(dustspecies, list):
            raise TypeError(f"{dustspecies} is not a list")

        for dspec in dustspecies:
            if dspec in cls._dust_species:
                logging.warning(f"{dspec} exists in the list of dust species, skip!")

            else:
                cls._dust_species.append(dspec)

        return cls._dust_species

    @classmethod
    def add_known_elements(cls, elements: list) -> list:
        """
        Add names of elements to the list of known elements

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
    def add_known_pseudoelements(cls, pelements: list) -> list:
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
                for e in Species._periodic_table + Species._isotopes_table
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
            basename = basename.replace(self._surface_prefix, "")

        # remove charge symbols
        if self.charge != 0:
            basename = re.sub(r"\+*$", "", basename)
            basename = re.sub(r"-*$", "", basename)

        return basename

    @property
    def binding_energy(self) -> float:
        """
        Binding energy of surface species.

        Return the binding energy if this species is at ice-phase. If the
        binding energy was not set before, it will be searched in the UMIST
        2012 binding energy data. Raise error when no value is found.

        Raises:
            RuntimeError: binding energy cannot be found
        """
        if self._binding_energy:
            return self._binding_energy

        if not self.is_surface:
            logging.fatal(f"{self.name} is not ice species, has no binding energy")

        self._binding_energy = (
            chemistry.user_binding_energy.get(self.name)
            if chemistry.user_binding_energy.get(self.name)
            else chemistry.rate12_binding_energy.get(self.gasname)
        )

        if not self._binding_energy:
            raise RuntimeError(f"Cannot find the binding energy of {self.name}")

        return self._binding_energy

    # alias function
    eb = binding_energy

    @property
    def charge(self) -> int:
        """
        Total charge of the species. Calculate the number of "+" and "-" in the end
        """
        if self.iselectron:
            return -1

        pcharge = "".join(re.findall(r"\+*$", self.name)).count("+")
        ncharge = "".join(re.findall(r"-*$", self.name)).count("-")
        return pcharge - ncharge

    @classmethod
    def dust_species(cls) -> list[str]:
        return cls._dust_species

    @property
    def enthalpy(self) -> float:
        """
        The enthalpy (in kJ/mol) of formation at 0K (surface species only, from CCCBDB)
        """
        if not self.is_surface:
            return 0.0

        if self._enthalpy:
            return self._enthalpy

        avails = {e.Species: float(e.Enthalpy) for e in Species._enthalpy_table}
        self._enthalpy = avails.get(self.gasname, 0.0)

        return self._enthalpy

    @property
    def gasname(self) -> str:
        """
        Return the name of its gas-phase species if this species is at ice-phase.
        Else return the current name.
        """
        return (
            self.name.replace(self._surface_prefix, "")
            if self.is_surface
            else self.name
        )

    @property
    def is_dust(self) -> bool:
        """
        Check whether the species is one form of dust

        Returns:
            bool: True if the species is one dust form
        """
        return self.name in self._dust_species

    @property
    def iselectron(self) -> bool:
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
        dust species (e.g. GRAIN0) if it is provided

        Returns:
            bool: True if the species is an atom
        """
        names, counts = self.element_count.keys(), self.element_count.values()
        return (
            len(names) == 1
            and sum(counts) == 1
            and self.charge == 0
            and not self.iselectron
            # and not self.is_dust
            and not self.is_surface
        )

    @property
    def is_surface(self) -> bool:
        """
        Check whether the species is sticking on surface (a surface species)
        """
        if self.is_dust:
            return False
        # if "GRAIN" in self.name.upper():
        #     return False
        return self.name.startswith(self._surface_prefix)

    # TODO: python 3.9 support classmethod property
    @classmethod
    def known_elements(cls) -> list:
        """
        Returns the current list of known elements

        Returns:
            list: the current list of known elements
        """
        return cls._known_elements

    @classmethod
    def known_pseudoelements(cls) -> list:
        """
        Returns the current list of known pseudo elements

        Returns:
            list: the current list of known pseudo elements
        """
        return cls._known_pseudoelements

    @property
    def mass(self) -> float:
        """
        The mass (amu) of the species, estimated by summing the mass of elements
        """
        if self._mass:
            return self._mass

        self._mass = 0.0
        for e in Species._periodic_table + Species._isotopes_table:
            self._mass += self._allcaps_element_count.get(e.Symbol.upper(), 0) * float(
                e.AtomicMass
            )

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
        for e in Species._periodic_table + Species._isotopes_table:
            self._massnumber += self._allcaps_element_count.get(e.Symbol.upper(), 0) * (
                float(e.NumberofNeutrons) + float(e.NumberofProtons)
            )

        # TODO: PAH, GRAIN, electron
        # self._massnumber += self.element_count.get("GRAIN", 0) * 1200.0

        if self._massnumber <= 0.0:
            logging.warning(f"{self.name} has zero or negative mass number!")

        return self._massnumber

    # alias function of mass number
    A = massnumber

    def photon_yield(self, default: float = 1e-3) -> float:
        """
        The photodesorption yield of the species (only for ice-phase species).
        Return default value if no value is found.

        Args:
            default (float, optional): default value of photodesorption yield.
                Defaults to 1e-3.

        Returns:
            float: photodesorption yield
        """
        if self._photon_yield:
            return self._photon_yield

        if not self.is_surface:
            logging.fatal(f"{self.name} is not ice species! No photodesorption yield")

        self._photon_yield = chemistry.user_photon_yield.get(self.name, default)

        return self._photon_yield

    @classmethod
    def remove_dust_species(cls, dustspecies: list[str]) -> list[str]:

        if not isinstance(dustspecies, list):
            raise TypeError(f"{dustspecies} is not a list")

        for dspec in dustspecies:
            cls._dust_species.remove(dspec)

        return cls._dust_species

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

        Reset known_elements, known_pseudoelements, dust_species, and surface_prefix.
        """
        cls._known_elements = []
        cls._known_pseudoelements = []
        cls._dust_species = []
        cls.surface_prefix = "#"

    @classmethod
    def set_dust_species(cls, dustspecies: list[str]) -> None:
        """Set the list of dust species

        Args:
            dustspecies (list[str]): the list of dust species
        """
        if not isinstance(dustspecies, list):
            raise TypeError(f"{dustspecies} is not a list")

        cls._dust_species = dustspecies.copy()

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
