import logging
import gac.singleton
from gac.auxiliary import request_user_check

default_element_list = [
    "e",
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

# default_element_list = [
#     "E",
#     "H",
#     "HE",
#     "LI",
#     "BE",
#     "B",
#     "C",
#     "N",
#     "O",
#     "F",
#     "NE",
#     "NA",
#     "MG",
#     "AL",
#     "SI",
#     "P",
#     "S",
#     "CL",
#     "AR",
#     "K",
#     "CA",
#     "SC",
#     "TI",
#     "V",
#     "CR",
#     "MN",
#     "FE",
#     "CO",
#     "NI",
#     "Cu",
#     "Zn",
#     "Ga",
#     "Ge",
# ]

default_pseudo_element_list = [
    "CR",
    "CRP",
    "Photon",
    "X",
    "p",
    "o",
    "c-",
    "l-",
]


def initialize_element(element_list):
    gac.singleton.element_list = element_list


def initialize_pseudo_element(pseudo_element_list):
    gac.singleton.pseudo_element_list = pseudo_element_list


def initialize_species(species_list):
    gac.singleton.species_list = species_list


def initialize_from_file(filename):
    with open(filename):
        pass


def initialize(filename=None):
    if not filename:
        logging.warning("No file assigned. Set Default species")
        initialize_element(default_element_list)
        initialize_pseudo_element(default_pseudo_element_list)
        # initialize_species(default_species_list)
    else:
        initialize_from_file(filename)

    all_element = gac.singleton.element_list + gac.singleton.pseudo_element_list

    dup_element = [ele for ele in all_element if all_element.count(ele) > 1]

    if len(dup_element) > 0:
        request_user_check(
            "{} are repeated in element list. Are you sure to continue? [Y/n]".format(
                dup_element
            ),
            "N",
        )
