import pytest
from pathlib import Path
from naunet.network import Network

inpath = Path("test/test_input")


@pytest.mark.slow
def test_init_network():
    network = Network()
    network.add_reaction_from_file(inpath / "deuspin.kida.uva.2017.in", "kida")
    network.check_duplicate_reaction()


def test_init_network_with_file():
    network = Network(inpath / "duplicate_test.dat", "kida")
    network = Network([inpath / "duplicate_test.dat"], "kida")
    network = Network(
        [inpath / "duplicate_test.dat", inpath / "simple_test.dat"], "kida"
    )
    network = Network(inpath / "duplicate_test.dat", ["kida"])
    network = Network([inpath / "duplicate_test.dat"], ["kida"])
    network = Network(
        [inpath / "duplicate_test.dat", inpath / "react_primordial.krome"],
        ["kida", "krome"],
    )


def test_check_duplicate():
    network = Network()
    network.add_reaction_from_file(inpath / "duplicate_test.dat", "kida")
    network.check_duplicate_reaction(full_check=False)


def test_generate_cvode_code_from_kida():
    network = Network()
    network.add_reaction_from_file(inpath / "simple_test.dat", "kida")
    network.check_duplicate_reaction()
    network.to_code(method="sparse", prefix="test/test_output/cvode_kida")


def test_generate_cvode_code_from_leeds():
    network = Network()
    network.add_reaction_from_file(inpath / "rate12_twobody_HO.rates", "leeds")
    network.check_duplicate_reaction()
    network.to_code(prefix="test/test_output/cvode_leeds")


def test_generate_cvode_code_from_krome():
    network = Network()
    network.add_reaction_from_file(inpath / "react_primordial.krome", "krome")
    network.to_code(prefix="test/test_output/cvode_krome")


def test_generate_odeint_code_from_krome():
    network = Network()
    network.add_reaction_from_file(inpath / "react_primordial.krome", "krome")
    prefix = "test/test_output/odeint_krome"
    network.to_code(solver="odeint", method="rosenbrock4", prefix=prefix)
