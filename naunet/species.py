import logging
import re
from collections import Counter
from . import settings


class Species:
    def __init__(self, name):
        self.name = name
        self.element_count = dict()
        self.charge = 0
        self.is_surface = False
        self._alias = None

        # initialize the default elements list if it is not
        # initialized when the fist species is instanciated
        if not settings.setting_initialized:
            settings.initialize()

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

    @staticmethod
    def _is_surface(name):
        return name[0] == settings.surface_symbol

    @staticmethod
    def _is_charged(name):
        return any([x in name[-1] for x in settings.charge_symbols])

    def _parse_molecule_name(self):
        element_sorted = sorted(
            settings.element_list + settings.pseudo_element_list, key=len, reverse=True
        )
        element_len = list(map(lambda x: len(x), element_sorted))

        specname = self.name
        # TODO: generalize the way to check surface and count charge
        if self._is_surface(specname):
            self.is_surface = True
            specname = specname.replace(settings.surface_symbol, "")

        if self._is_charged(specname):
            pcharge = "".join(re.findall(r"\+*$", specname)).count("+")
            ncharge = "".join(re.findall(r"-*$", specname)).count("-")
            self.charge = pcharge - ncharge
            # remove the charge symbols
            specname = re.sub(r"\+*$", "", specname)
            specname = re.sub(r"-*$", "", specname)

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
    def basename(self) -> str:
        """
        Clip the surface and charge symbols. Return the name of the molecular/atom

        Returns:
            str: name of the molecular/atom without charge or phase information
        """
        # TODO: generalize the way to check surface and count charge
        basename = self.name
        if self.is_surface:
            basename = basename.replace(settings.surface_symbol, "")

        if self.charge > 0:
            basename = re.sub(r"\+*$", "", basename)
        elif self.charge < 0:
            basename = re.sub(r"-*$", "", basename)

        return basename

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
