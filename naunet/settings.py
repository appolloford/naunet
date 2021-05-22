"""Global settings"""

import logging
from .auxiliary import conti_confirm

setting_initialized = False

element_list = []

pseudo_element_list = []

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

default_pseudo_element_list = [
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


def add_element(elements: list) -> list:
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
    return element_list


def add_pseudo_element(elements: list) -> list:
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
    return pseudo_element_list


def duplicate_elements() -> list:
    all_element = element_list + pseudo_element_list
    dup_element = [ele for ele in all_element if all_element.count(ele) > 1]
    return dup_element


def _initialize_element(input_element_list: list) -> None:
    element_list.clear()
    element_list.extend(input_element_list)


def _initialize_pseudo_element(input_pseudo_element_list: list) -> None:
    pseudo_element_list.clear()
    pseudo_element_list.extend(input_pseudo_element_list)


def initialize(
    element: list = None,
    pseudo_element: list = None,
):

    global setting_initialized

    if element or pseudo_element:

        _initialize_element(element)
        _initialize_pseudo_element(pseudo_element)

    else:

        logging.warning("No assigned data. Use default elements")

        _initialize_element(default_element_list)
        _initialize_pseudo_element(default_pseudo_element_list)

    dup_element = duplicate_elements()

    if dup_element:
        conti_confirm(
            "{} are repeated in element list. Continue?".format(dup_element),
            default=True,
        )

    setting_initialized = True


def remove_element(elements: list) -> list:
    for ele in elements:
        element_list.remove(ele)
    return element_list


def remove_pseudo_element(elements: list) -> list:
    for ele in elements:
        pseudo_element_list.remove(ele)
    return pseudo_element_list


def reset_element_list(new_element_list: list = None) -> list:
    element_list.clear()
    element_list.extend(new_element_list)
    return element_list


def reset_pseudo_element_list(new_pseudo_element_list: list = None) -> list:
    pseudo_element_list.clear()
    pseudo_element_list.extend(new_pseudo_element_list)
    return pseudo_element_list
