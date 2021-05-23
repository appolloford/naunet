import logging
import re
import csv
import os
from collections import namedtuple
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
        "g",
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
    ]
    surface_prefix = "#"

    _known_elements = []
    _known_pseudoelements = []
    _periodic_table = []

    def __init__(self, name):
        self.name = name
        self.element_count = dict()
        self._alias = None
        self._binding_energy = None
        self._photon_yield = None
        self._mass = 0.0
        self._massnumber = 0

        # Initialize known elements if not set when the fist species is instanciated
        if not Species._known_elements and not Species._known_pseudoelements:
            logging.warning("No assigned element list. Use default elements")

            Species._known_elements.extend(Species.default_elements)
            Species._known_pseudoelements.extend(Species.default_pseudoelements)

        # Create periodic table if not created
        if not Species._periodic_table:
            path = os.path.dirname(chemistry.__file__)
            # Load the periodic table as a list of namedtuple
            # Source: https://gist.github.com/GoodmanSciences/c2dd862cd38f21b0ad36b8f96b4bf1ee
            # TODO: isotopes
            with open(os.path.join(path, "periodictable.csv"), newline="") as infile:
                # remove comment lines in case
                reader = csv.reader(filter(lambda row: row[0] != "#", infile))
                Element = namedtuple("Element", next(reader))
                Species._periodic_table = list(map(Element._make, reader))

        self._parse_molecule_name()

    def __eq__(self, o: object) -> bool:
        if isinstance(o, Species):
            return self.name == o.name
        return NotImplemented

    def __hash__(self) -> int:
        return hash(self.name)

    def __lt__(self, o) -> bool:
        if isinstance(o, Species):
            return self.name < o.name
        return NotImplemented

    def __str__(self) -> str:
        return "Species({})".format(self.name)

    def _add_element_count(self, element, count):
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
        all_element = cls._known_elements + cls._known_pseudoelements
        dup_element = [ele for ele in all_element if all_element.count(ele) > 1]
        if dup_element:
            logging.error(
                f"Repeated definitions of {dup_element} in elements and pseudoelemts"
            )

    def _parse_molecule_name(self):
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
    def add_known_elements(cls, elements: list) -> list:
        for ele in elements:
            if ele in cls._known_elements:
                logging.warning("{} exists in element list, skip!".format(ele))
            elif ele in cls._known_pseudoelements:
                logging.warning(
                    "{} exists in pseudo element list, move to element list!".format(
                        ele
                    )
                )
                cls._known_pseudoelements.remove(ele)
                cls._known_elements.append(ele)
            else:
                cls._known_elements.append(ele)
        cls._check_elements()
        return cls._known_elements

    @classmethod
    def add_known_pseudoelements(cls, pelements: list) -> list:
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
        str: Alias of the species.

        Replace surface symbol by `G` and charges by `M`(-), `I`(+). Used in
        indexing species in codes. Can be customized by setter if needed.
        """

        if not self._alias:
            basename = self.basename
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
        """str: The species name without surface symbols and charges."""

        basename = self.name
        # remove surface symbol
        if self.is_surface:
            basename = basename.replace(Species.surface_prefix, "")

        # remove charge symbols
        if self.charge != 0:
            basename = re.sub(r"\+*$", "", basename)
            basename = re.sub(r"-*$", "", basename)

        return basename

    @property
    def binding_energy(self) -> float:
        """
        float: Binding energy of surface species.

        Return the binding energy if this species is at ice-phase. If the
        binding energy was not set before, it will be searched in the UMIST
        2012 binding energy data. Raise error when no value is found.

        Raises:
            RuntimeError: binding energy cannot be found in UMIST2012 or
            customized data
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

    @property
    def charge(self) -> int:
        """
        Total charge of the species. Calculate the number of "+" and "-" in the end

        Returns:
            int: total charge
        """
        pcharge = "".join(re.findall(r"\+*$", self.name)).count("+")
        ncharge = "".join(re.findall(r"-*$", self.name)).count("-")
        return pcharge - ncharge

    @property
    def gasname(self) -> str:
        """
        Return the name of its gas-phase species if this species is at ice-phase.
        Else return the current name.

        Returns:
            str: name of its gas-phase species
        """
        return (
            self.name.replace(Species.surface_prefix, "")
            if self.is_surface
            else self.name
        )

    @property
    def is_surface(self):
        if "GRAIN" in self.name.upper():
            return False
        return self.name.startswith(Species.surface_prefix)

    # TODO: python 3.9 support classmethod property
    @classmethod
    def known_elements(cls) -> list:
        return cls._known_elements

    @classmethod
    def known_pseudoelements(cls) -> list:
        return cls._known_pseudoelements

    @property
    def mass(self) -> float:
        """
        The mass (amu) of the species, estimated by summing the mass of elements

        Returns:
            float: mass of the species
        """
        if self._mass:
            return self._mass

        for e in Species._periodic_table:
            self._mass += self.element_count.get(e.Symbol, 0) * float(e.AtomicMass)

        # ? electron mass
        # self._mass -= 0.00054858 * self.charge

        # TODO: GRAIN, electron and isotopes

        return self._mass

    @property
    def massnumber(self) -> float:
        """
        The mass number (neutron + proton) of the species

        Returns:
            int: mass number of the species
        """
        if self._massnumber:
            return self._massnumber

        for e in Species._periodic_table:
            self._massnumber += self.element_count.get(e.Symbol, 0) * (
                float(e.NumberofNeutrons) + float(e.NumberofProtons)
            )

        # TODO: GRAIN, electron and isotopes
        self._massnumber += self.element_count.get("D", 0) * 2.0
        self._massnumber += self.element_count.get("GRAIN", 0) * 1200.0

        return self._massnumber

    @property
    def photon_yield(self) -> float:
        """
        The photodesorption yield of the species (ice-phase only). Return default value
        (1.0e-3) if no value is found. The default value can be changed in
        naunet.chemistry.default_photon_yield

        Returns:
            float: photodesorption yield
        """
        if self._photon_yield:
            return self._photon_yield

        if not self.is_surface:
            logging.fatal(f"{self.name} is not ice species! No photodesorption yield")

        self._photon_yield = chemistry.user_photon_yield.get(
            self.name, chemistry.default_photon_yield
        )

        return self._photon_yield

    @classmethod
    def remove_known_elements(cls, elements: list) -> list:
        for ele in elements:
            cls._known_elements.remove(ele)
        cls._check_elements()
        return cls._known_elements

    @classmethod
    def remove_known_pseudoelements(cls, pelemets: list) -> list:
        for ele in pelemets:
            cls._known_pseudoelements.remove(ele)
        cls._check_elements()
        return cls._known_pseudoelements

    @classmethod
    def reset(cls) -> None:
        """Reset class attributes in Species

        Reset known_elements, known_pseudoelements, and surface_prefix.
        """
        cls._known_elements = []
        cls._known_pseudoelements = []
        cls.surface_prefix = "#"

    @classmethod
    def set_known_elements(cls, elements: list) -> None:
        cls._known_elements.clear()
        cls._known_elements.extend(elements)
        cls._check_elements()

    @classmethod
    def set_known_pseudoelements(cls, pelements: list) -> None:
        cls._known_pseudoelements.clear()
        cls._known_pseudoelements.extend(pelements)
        cls._check_elements()


def top_abundant_species(species_list, abundances, element=None, rank=-1):
    """
    The function returns a tuple list sorted by the abundances of element.

    Args:
        species_list (list): list of `Molecule()` objects.
        abundances (list): The abundances of the species in the species_list.
        element (string, optional): The target element. The order is sorted by the abundances weighted by the number of element in the species. Defaults to None.
        rank (int, optional): Return the species in the top number. Defaults to -1 (all sorted species).

    Raises:
        RuntimeError: The element could doesn't exist in the species_list and return an empty list.

    Returns:
        list: Tuple of `(Molecule, float)`.
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
