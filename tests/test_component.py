from typing import OrderedDict
from naunet.component import Component, VariableType as vt


class Child(Component):
    def __init__(self) -> None:
        super().__init__()

        self.register("density", ("nH", None, vt.param))
        self.register("eb_h2d", ("eb_h2d", 1.21e3, vt.constant))
        self.register("grain_density", ("gdens", "y[GRAIN0]", vt.derived))


def test_child():
    child = Child()
    child.symbols.density.symbol == "nH"
    child.symbols.density.value == "nH"
    child.symbols.density.type == vt.param
    child.constants == {"eb_h2d": 1.21e3}
    deriveds = OrderedDict()
    deriveds["gdens"] = "y[GRAIN0]"
    child.deriveds == deriveds
