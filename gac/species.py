import logging
from math import sin
import re
from collections import Counter
from gac import singleton
from gac.singleton import (
    element_list,
    pseudo_element_list,
    surface_symbol,
    charge_symbols,
)

# modify the local copy to contain the pseudo elements
element_list = element_list + pseudo_element_list


class Species:
    def __init__(self, name):
        self.name = name
        self.element_count = dict()
        self.is_surface = False

        self._parse_molecule_name()

    def __eq__(self, o: object) -> bool:
        if isinstance(o, Species):
            return self.name == o.name
        return NotImplemented

    def __str__(self) -> str:
        return "Species({})".format(self.name)

    def _add_element_count(self, element, count):
        if element in pseudo_element_list:
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
        return name[0] == surface_symbol

    @staticmethod
    def _is_charged(name):
        return any([x in name[-1] for x in charge_symbols])

    def _parse_molecule_name(self):
        element_sorted = sorted(element_list, key=len, reverse=True)
        element_len = list(map(lambda x: len(x), element_sorted))

        specname = self.name
        # TODO: generalize the way to check surface and count charge
        if self._is_surface(specname):
            self.is_surface = True
            specname = specname.replace(surface_symbol, "")

        if self._is_charged(specname):
            pcharge = "".join(re.findall(".+$", specname)).count("+")
            ncharge = "".join(re.findall(".-$", specname)).count("-")
            self.charge = pcharge - ncharge
            # remove the charge symbols
            specname = re.sub(".+$", "", specname)
            specname = re.sub(".-$", "", specname)

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
