from types import new_class
from naunet import network
from naunet.network import Network


def test_init_network():
    network = Network("test/deuspin.kida.uva.2017.in", "kida")
    network.check_duplicate_reaction()


def test_check_duplicate():
    network = Network("test/duplicate_test.dat", "kida")
    network.check_duplicate_reaction(full_check=False)


def test_generate_ccode():
    network = Network("test/simple_test.dat", "kida")
    network.check_duplicate_reaction()
    network.collect_infos()
    network.constants_header()
    for func in ["fex", "jtv", "jac"]:
        network.create_ode_expression().to_ccode(
            function=func,
            use_template=True,
            to_file=True,
            file_name=func,
            prefix="./",
            header=True,
        )
