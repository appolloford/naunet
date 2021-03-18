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
    network.info.to_ccode(
        to_file=True,
        prefix="test/test_output/",
        file_name="naunet_constants.h",
    )
    network.ode_expression.to_ccode(
        to_file=True,
        prefix="test/test_output/",
        file_name=f"naunet_ode.cpp",
        header=True,
        header_prefix="test/test_output/",
        header_file=f"naunet_ode.h",
    )


def test_generate_ccode_from_leeds():
    network = Network("test/test_input/rate12_twobody_HO.rates", "leeds")
    network.check_duplicate_reaction()
    network.info.to_ccode(
        to_file=True,
        prefix="test/test_output/",
        file_name="naunet_constants.h",
    )
    network.ode_expression.to_ccode(
        to_file=True,
        prefix="test/test_output/",
        file_name=f"naunet_ode.cpp",
        header=True,
        header_prefix="test/test_output/",
        header_file=f"naunet_ode.h",
    )


def test_generate_ccode_from_krome():
    network = Network("test/test_input/react_primordial.krome", "krome")
    network.info.to_ccode(
        to_file=True,
        prefix="test/test_output/",
        file_name="naunet_constants.h",
    )
    network.userdata.to_ccode(
        to_file=True,
        prefix="test/test_output/",
        file_name="naunet_userdata.h",
    )
    network.ode_expression.to_ccode(
        to_file=True,
        prefix="test/test_output/",
        file_name=f"naunet_ode.cpp",
        header=True,
        header_prefix="test/test_output/",
        header_file=f"naunet_ode.h",
    )