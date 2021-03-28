from sys import prefix
import pytest
from naunet.network import Network


@pytest.mark.slow
def test_init_network():
    network = Network("test/test_input/deuspin.kida.uva.2017.in", "kida")
    network.check_duplicate_reaction()


def test_check_duplicate():
    network = Network("test/test_input/duplicate_test.dat", "kida")
    network.check_duplicate_reaction(full_check=False)


def test_generate_ccode_from_kida():
    network = Network("test/test_input/simple_test.dat", "kida")
    network.check_duplicate_reaction()
    network.to_code(linsolver="sparse", prefix="test/test_output/")
    network.finalize()


def test_generate_ccode_from_leeds():
    network = Network("test/test_input/rate12_twobody_HO.rates", "leeds")
    network.check_duplicate_reaction()
    network.to_code(prefix="test/test_output/")
    network.finalize()


def test_generate_ccode_from_krome():
    network = Network("test/test_input/react_primordial.krome", "krome")
    network.to_code(prefix="test/test_output/")
    network.finalize()