import logging
import re
import csv
import os
from collections import namedtuple
from . import settings
from . import chemistry

element_initialized = False
element_table = []


def initialize_elements():
    path = os.path.dirname(settings.__file__)
    global element_initialized, element_table
    # Load the periodic table as a list of namedtuple
    # Source: https://gist.github.com/GoodmanSciences/c2dd862cd38f21b0ad36b8f96b4bf1ee
    # TODO: isotopes
    with open(os.path.join(path, "chemistry/periodictable.csv"), newline="") as infile:
        # remove comment lines in case
        reader = csv.reader(filter(lambda row: row[0] != "#", infile))
        Element = namedtuple("Element", next(reader))
        element_table = list(map(Element._make, reader))
    element_initialized = True


class Species:
    def __init__(self, name):
        self.name = name
        self.element_count = dict()
        self._alias = None
        self._binding_energy = None
        self._photon_yield = None
        self._mass = 0.0
        self._massnumber = 0

        # initialize the default elements list if it is not
        # initialized when the fist species is instanciated
        if not settings.setting_initialized:
            settings.initialize()

        # Create table of elements if not created
        if not element_initialized:
            initialize_elements()

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
        if element in settings.pseudo_element_list:
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

    def _parse_molecule_name(self):
        element_sorted = sorted(
            settings.element_list + settings.pseudo_element_list, key=len, reverse=True
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

    @property
    def alias(self) -> str:
        """
        Alias of the species. Special symbols are removed. Used in indexing species in codes.

        Returns:
            str: alias of the species
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
    def alias(self, name):
        """
        Setter of alias. Customized the alias if needed.

        Args:
            name (str): new alias of the species
        """
        self._alias = name

    @property
    def basename(self) -> str:
        """
        Clip the surface and charge symbols. Return the name of the molecular/atom

        Returns:
            str: name of the molecular/atom without charge or phase information
        """
        basename = self.name

        # remove surface symbol
        if self.is_surface:
            basename = basename.replace(settings.surface_symbol, "")

        # remove charge symbols
        if self.charge != 0:
            basename = re.sub(r"\+*$", "", basename)
            basename = re.sub(r"-*$", "", basename)

        return basename

    @property
    def binding_energy(self) -> float:
        """
        Return the binding energy if this species is at ice-phase. If the
        binding energy was not set before, it will be searched in the UMIST
        2012 binding energy data. Raise error when no value is found.

        Returns:
            float: the binding energy
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
            self.name.replace(settings.surface_symbol, "")
            if self.is_surface
            else self.name
        )

    @property
    def is_surface(self):
        if "GRAIN" in self.name.upper():
            return False
        return self.name.startswith(settings.surface_symbol)

    @property
    def mass(self) -> float:
        """
        The mass (amu) of the species, estimated by summing the mass of elements

        Returns:
            float: mass of the species
        """
        if self._mass:
            return self._mass

        for e in element_table:
            self._mass += self.element_count.get(e.Symbol, 0) * float(e.AtomicMass)

        # ? electron mass
        # self._mass -= 0.00054858 * self.charge

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

        for e in element_table:
            self._massnumber += self.element_count.get(e.Symbol, 0) * (
                float(e.NumberofNeutrons) + float(e.NumberofProtons)
            )

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
