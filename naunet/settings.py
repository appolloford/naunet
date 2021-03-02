"""
This module is provided as a singleton. Modify naunet.settings.xxx to share the values globally
"""
import sympy as sym
import logging
from .auxiliary import request_user_check

setting_initialized = False

species_initialized = False

element_list = []

pseudo_element_list = []

species_list = []

surface_symbol = "#"

charge_symbols = ["+", "-"]

user_symbols = {}

default_user_symbols = {
    "CRIR": sym.symbols("zeta"),
    "Temperature": sym.symbols("Tgas"),
    "VisualExtinction": sym.symbols("Av"),
    "UVPHOT": sym.symbols("uv"),
}

ode_symbols = {
    "ode_vector": "y",
    "fex_lhs": "ydot",
    "jtv_lhs": "jv",
    "jac_lhs": "",  # dummy, LHS of jac() is hard coded
    "rate": "k",
}

default_element_list = [
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
    "XRAY",
    "Photon",
    "X",
    "M",
    "p",
    "o",
    "m",
    "c-",
    "l-",
    "*",
]


def add_element(elements: list):
    for ele in elements:
        if ele in element_list:
            logging.warning("{} exists in element list, skip!".format(ele))
        elif ele in pseudo_element_list:
            logging.warning(
                "{} exists in pseudo element list, move to element list!".format(ele)
            )
            pseudo_element_list.remove(ele)
            element_list.append(ele)
        else:
            element_list.append(ele)


def add_pseudo_element(elements: list):
    for ele in elements:
        if ele in pseudo_element_list:
            logging.warning("{} exists in pseudo element list, skip!".format(ele))
        elif ele in element_list:
            logging.warning(
                "{} exists in element list, move to pseudo element list!".format(ele)
            )
            element_list.remove(ele)
            pseudo_element_list.append(ele)
        else:
            pseudo_element_list.append(ele)


def _initialize_element(input_element_list):
    element_list.clear()
    element_list.extend(input_element_list)


def _initialize_pseudo_element(input_pseudo_element_list):
    pseudo_element_list.clear()
    pseudo_element_list.extend(input_pseudo_element_list)


def _initialize_species(input_species_list):
    species_list.clear()
    species_list.extend(input_species_list)


def _initialize_user_symbols(input_user_symbols):
    user_symbols.update(input_user_symbols)


def _initialize_from_file(filename):
    with open(filename):
        pass


def initialize(
    element: list = [],
    pseudo_element: list = [],
    species: list = [],
    filename: str = None,
):

    global species_initialized, setting_initialized

    if filename:

        _initialize_from_file(filename)

    elif len(element) or len(pseudo_element):

        _initialize_element(element)
        _initialize_pseudo_element(pseudo_element)

        if len(species):
            _initialize_species(species)
            species_initialized = True

    else:

        logging.warning("No assigned data. Use default elements")

        # print(default_pseudo_element_list, default_element_list)
        _initialize_element(default_element_list)
        _initialize_pseudo_element(default_pseudo_element_list)
        _initialize_user_symbols(default_user_symbols)
        # initialize_species(default_species_list)

    all_element = element_list + pseudo_element_list

    dup_element = [ele for ele in all_element if all_element.count(ele) > 1]

    if len(dup_element) > 0:
        request_user_check(
            "{} are repeated in element list. Are you sure to continue? [Y/n]".format(
                dup_element
            ),
            "N",
        )

    user_symbols.update(default_user_symbols)

    setting_initialized = True


def remove_element(elements):
    for ele in elements:
        element_list.remove(ele)


def remove_pseudo_element(elements):
    for ele in elements:
        pseudo_element_list.remove(ele)


def reset_element_list(new_element_list: list = []):
    element_list.clear()
    element_list.extend(new_element_list)


def reset_pseudo_element_list(new_pseudo_element_list: list = []):
    pseudo_element_list.clear()
    pseudo_element_list.extend(new_pseudo_element_list)


def check_duplicate_species():
    all_element = element_list + pseudo_element_list
    dup_element = [ele for ele in all_element if all_element.count(ele) > 1]
    if len(dup_element) > 0:
        print("These elements are duplicate: ", dup_element)
