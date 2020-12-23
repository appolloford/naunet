from gac.settings import (
    element_list,
    default_element_list,
    add_element,
    reset_element_list,
)


def test_module():
    assert element_list == default_element_list
    add_element(["Ne"])
    assert element_list == default_element_list + ["Ne"]
