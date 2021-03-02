from naunet.settings import (
    element_list,
    default_element_list,
    add_element,
    reset_element_list,
    initialize,
)


def test_module():
    initialize()
    assert element_list == default_element_list
    add_element(["Ne"])
    assert element_list == default_element_list + ["Ne"]
