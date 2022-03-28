import pytest
from pathlib import Path
from naunet.species import Species
from naunet.reactions.reaction import Reaction
from naunet.reactions.kidareaction import KIDAReaction
from naunet.network import Network

inpath = Path("test/test_input")


def test_init_network():
    network = Network()

    assert network.format_list == set()
    assert network.reaction_list == []
    assert network.reactants_in_network == set()
    assert network.products_in_network == set()


def test_add_reaction():
    reac = Reaction(
        ["He", "CR"],
        ["He+", "e-"],
        -9999.0,
        9999.0,
        0.5,
        0.0,
        0.0,
    )
    network = Network()

    network.add_reaction(reac)

    assert network.format_list == set()
    assert network.reaction_list == [reac]
    assert network.reactants_in_network == {Species("He")}
    assert network.products_in_network == {Species("He+"), Species("e-")}

    react_string = "".join(
        [
            "He         CR                     ",
            "He+        e-                                            ",
            "5.000e-01  0.000e+00  0.000e+00 ",
            "2.00e+00 0.00e+00 ",
            "logn  1  -9999   9999  1     3 1  1 ",
        ]
    )

    network = Network()
    network.add_reaction((react_string, "kida"))

    assert network.format_list == set(["kida"])
    assert network.reaction_list == [KIDAReaction(react_string)]
    assert network.reactants_in_network == {Species("He")}
    assert network.products_in_network == {Species("He+"), Species("e-")}


def test_init_network_from_reactions():
    reactions = []
    reactions.append(
        Reaction(
            ["H", "CR"],
            ["H+", "e-"],
            10.0,
            41000.0,
            5.98e-18,
            0.0,
            0.0,
        )
    )

    reactions.append(
        Reaction(
            ["H2", "CH"],
            ["C", "H2", "H"],
            1340.0,
            41000.0,
            6e-9,
            0.0,
            40200.0,
            format="kida",
        )
    )

    reactions.append(
        Reaction(
            ["H2", "e-"],
            ["H", "H", "e-"],
            3400.0,
            41000.0,
            3.22e-09,
            0.35,
            102000.0,
            format="umist",
        )
    )

    network = Network(reactions)
    assert network.format_list == set(["kida", "umist"])
    assert network.reaction_list == reactions
    assert network.reactants_in_network == {
        Species("H"),
        Species("H2"),
        Species("CH"),
        Species("e-"),
    }
    assert network.products_in_network == {
        Species("H"),
        Species("H+"),
        Species("H2"),
        Species("C"),
        Species("e-"),
    }


def test_init_network_with_file():
    network = Network(
        filelist=inpath / "duplicate_test.dat",
        fileformats="kida",
    )
    network = Network(
        filelist=[inpath / "duplicate_test.dat"],
        fileformats="kida",
    )
    network = Network(
        filelist=[inpath / "duplicate_test.dat", inpath / "simple_test.dat"],
        fileformats="kida",
    )
    network = Network(
        filelist=inpath / "duplicate_test.dat",
        fileformats=["kida"],
    )
    network = Network(
        filelist=[inpath / "duplicate_test.dat"],
        fileformats=["kida"],
    )
    network = Network(
        filelist=[inpath / "duplicate_test.dat", inpath / "react_primordial.krome"],
        fileformats=["kida", "krome"],
    )


def test_find_species():

    reactions = []
    reactions.append(
        Reaction(
            ["H", "CR"],
            ["H+", "e-"],
            10.0,
            41000.0,
            5.98e-18,
            0.0,
            0.0,
        )
    )

    reactions.append(
        Reaction(
            ["H2", "CH"],
            ["C", "H2", "H"],
            1340.0,
            41000.0,
            6e-9,
            0.0,
            40200.0,
        )
    )

    reactions.append(
        Reaction(
            ["H2", "e-"],
            ["H", "H", "e-"],
            3400.0,
            41000.0,
            3.22e-09,
            0.35,
            102000.0,
        )
    )

    network = Network(reactions)
    assert network.find_reactant("H2") == [1, 2]
    assert network.find_product("e-") == [0, 2]
    assert network.find_species("C") == [1]
    assert network.find_species("H") == [0, 1, 2]


@pytest.mark.slow
def test_init_network_from_kida():
    network = Network()
    network.add_reaction_from_file(inpath / "deuspin.kida.uva.2017.in", "kida")
    network.check_duplicate_reaction()


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


def test_write_network():
    network = Network(
        filelist=inpath / "duplicate_test.dat",
        fileformats="kida",
    )
    network.write("test/test_output/writereaction.txt")
