from types import new_class
from gac import network
from gac.network import Network


def test_init_network():
    # network = Network("test/deuspin.kida.uva.2017.in", "kida")
    network = Network("test/simple_test.dat", "kida")
    network.check_duplicate_reaction()
    network.to_ccode()