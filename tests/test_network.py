import os
import subprocess
import pytest
from pathlib import Path
from naunet.species import Species
from naunet.reactions.reaction import Reaction
from naunet.reactions.reactiontype import ReactionType
from naunet.reactions.kidareaction import KIDAReaction
from naunet.network import Network

GITHUB_ACTIONS = os.getenv("GITHUB_ACTIONS") == "true"


@pytest.fixture
def example_reaction1():
    return Reaction(
        ["He", "CR"],
        ["He+", "e-"],
        -9999.0,
        9999.0,
        0.5,
        0.0,
        0.0,
        reaction_type=ReactionType.GAS_COSMICRAY,
    )


@pytest.fixture
def example_reaction2():
    return Reaction(
        ["N2", "Photon"],
        ["N", "N"],
        -9999.0,
        9999.0,
        2.3e-10,
        0.0,
        3.88,
        reaction_type=ReactionType.GAS_PHOTON,
    )


@pytest.fixture
def example_reaction_string1():
    return "".join(
        [
            "He         CR                     ",
            "He+        e-                                            ",
            "5.000e-01  0.000e+00  0.000e+00 ",
            "2.00e+00 0.00e+00 ",
            "logn  1  -9999   9999  1     3 1  1 ",
        ]
    )


@pytest.fixture
def example_reaction_string2():
    return "".join(
        [
            "N2         Photon                 ",
            "N          N                                             ",
            "2.300e-10  0.000e+00  3.880e+00 ",
            "2.00e+00 0.00e+00 ",
            "logn  3  -9999   9999  2   364 1  1",
        ]
    )


@pytest.fixture
def example_reaction_list():
    return [
        Reaction(
            ["H", "CR"],
            ["H+", "e-"],
            10.0,
            41000.0,
            5.98e-18,
            0.0,
            0.0,
            ReactionType.GAS_COSMICRAY,
        ),
        Reaction(
            ["H2", "CH"],
            ["C", "H2", "H"],
            1340.0,
            41000.0,
            6e-9,
            0.0,
            40200.0,
            ReactionType.GAS_TWOBODY,
            format="kida",
        ),
        Reaction(
            ["H2", "e-"],
            ["H", "H", "e-"],
            3400.0,
            41000.0,
            3.22e-09,
            0.35,
            102000.0,
            ReactionType.GAS_TWOBODY,
            format="umist",
        ),
    ]


@pytest.fixture
def empty_network():
    return Network()


@pytest.fixture
def example_network_from_reaction_list(example_reaction_list):
    return Network(example_reaction_list)


def test_init_network(empty_network):
    assert empty_network.format_list == set()
    assert empty_network.reaction_list == []
    assert empty_network._reactants == set()
    assert empty_network._products == set()


def test_add_reaction(empty_network, example_reaction1):
    empty_network.add_reaction(example_reaction1)

    assert empty_network.format_list == {"naunet"}
    assert empty_network.reaction_list == [example_reaction1]
    assert empty_network._reactants == {Species("He")}
    assert empty_network._products == {Species("He+"), Species("e-")}


def test_add_reaction_string(empty_network, example_reaction_string1):
    empty_network.add_reaction((example_reaction_string1, "kida"))
    assert empty_network.format_list == set(["kida"])
    assert empty_network.reaction_list == [KIDAReaction(example_reaction_string1)]
    assert empty_network._reactants == {Species("He")}
    assert empty_network._products == {Species("He+"), Species("e-")}


def test_init_network_from_reactions(
    example_network_from_reaction_list,
    example_reaction_list,
):
    network = example_network_from_reaction_list
    assert network.format_list == set(["kida", "umist", "naunet"])
    assert network.reaction_list == example_reaction_list
    assert network._reactants == {
        Species("H"),
        Species("H2"),
        Species("CH"),
        Species("e-"),
    }
    assert network._products == {
        Species("H"),
        Species("H+"),
        Species("H2"),
        Species("C"),
        Species("e-"),
    }


def test_init_network_with_file(datadir):
    network = Network(
        filelist=datadir / "duplicate.kida",
        fileformats="kida",
    )
    network = Network(
        filelist=[datadir / "duplicate.kida"],
        fileformats="kida",
    )
    network = Network(
        filelist=[datadir / "duplicate.kida", datadir / "minimal.kida"],
        fileformats="kida",
    )
    network = Network(
        filelist=datadir / "duplicate.kida",
        fileformats=["kida"],
    )
    network = Network(
        filelist=[datadir / "duplicate.kida"],
        fileformats=["kida"],
    )
    network = Network(
        filelist=[datadir / "duplicate.kida", datadir / "primordial.krome"],
        fileformats=["kida", "krome"],
    )


def test_find_species(
    example_network_from_reaction_list,
    example_reaction_list,
):

    network = example_network_from_reaction_list
    assert network.where(example_reaction_list[0]) == [0]
    assert network.where(species="H2", mode="reactant") == [1, 2]
    assert network.where(species="e-", mode="product") == [0, 2]
    assert network.where(species="C") == [1]
    assert network.where(species="H") == [0, 1, 2]


@pytest.mark.slow
def test_init_network_from_kida(datadir):
    network = Network()
    network.add_reaction_from_file(datadir / "deuspin.kida", "kida")
    network.check_duplicate_reaction()


def test_check_duplicate(datadir):
    network = Network()
    network.add_reaction_from_file(datadir / "duplicate.kida", "kida")
    network.check_duplicate_reaction(full_check=False)


def test_generate_cvode_code_from_kida(tmp_path, datadir):
    network = Network()
    network.add_reaction_from_file(datadir / "minimal.kida", "kida")
    network.check_duplicate_reaction()
    network.to_code(method="sparse", prefix=tmp_path / "cvode_kida")


def test_generate_cvode_code_from_leeds(tmp_path, datadir):
    network = Network()
    network.add_reaction_from_file(datadir / "rate12_HO.leeds", "leeds")
    network.check_duplicate_reaction()
    network.to_code(prefix=tmp_path / "cvode_leeds")


def test_generate_cvode_code_from_krome(tmp_path, datadir):
    network = Network()
    network.add_reaction_from_file(datadir / "primordial.krome", "krome")
    network.to_code(prefix=tmp_path / "cvode_krome")


def test_generate_odeint_code_from_krome(tmp_path, datadir):
    network = Network()
    network.add_reaction_from_file(datadir / "primordial.krome", "krome")
    prefix = tmp_path / "odeint_krome"
    network.to_code(solver="odeint", method="rosenbrock4", prefix=prefix)


def test_write_network(tmp_path, datadir):
    network = Network(
        filelist=datadir / "duplicate.kida",
        fileformats="kida",
    )
    network.write(tmp_path / "reactions.txt")
    network.write(tmp_path / "reactions.naunet", format="naunet")
    network.write(tmp_path / "reactions.krome", format="krome")


@pytest.mark.skipif(GITHUB_ACTIONS, reason="on github")
def test_export_empty_network(tmp_path):
    network = Network()
    network.export(prefix=tmp_path / "empty_network", overwrite=True)

    os.chdir(tmp_path / "empty_network")
    process = subprocess.Popen(
        ["cmake", "-S", ".", "-B", "build"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    stdout, stderr = process.communicate()

    assert "CMake Error" not in stderr.decode("utf-8")

    # test the export code can be compiled
    process = subprocess.Popen(
        ["cmake", "--build", "build"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    stdout, stderr = process.communicate()

    assert " failed" not in stdout.decode("utf-8")
    assert " error:" not in stderr.decode("utf-8")

    os.chdir("../../../")


# TODO: update the github workflow to have proper environment for the test
@pytest.mark.skipif(GITHUB_ACTIONS, reason="on github")
def test_export_network(tmp_path, datadir):
    network = Network(
        filelist=datadir / "minimal.kida",
        fileformats="kida",
    )
    # network = Network(
    #     filelist=datadir / "primordial.krome",
    #     fileformats="krome",
    # )
    network.export(prefix=tmp_path / "network_export", overwrite=True)

    os.chdir(tmp_path / "network_export")
    process = subprocess.Popen(
        ["cmake", "-S", ".", "-B", "build"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    stdout, stderr = process.communicate()

    assert "CMake Error" not in stderr.decode("utf-8")

    # test the export code can be compiled
    process = subprocess.Popen(
        ["cmake", "--build", "build"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    stdout, stderr = process.communicate()

    assert " failed" not in stdout.decode("utf-8")
    assert " error:" not in stderr.decode("utf-8")

    # check the config file is valid
    process = subprocess.Popen(
        ["naunet", "render", "-f"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    stdout, stderr = process.communicate()

    assert "Error" not in stdout.decode("utf-8")

    # check the re-renderred result still can be compiled
    process = subprocess.Popen(
        ["cmake", "--build", "build"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    stdout, stderr = process.communicate()

    assert " failed" not in stdout.decode("utf-8")
    assert " error:" not in stderr.decode("utf-8")

    # os.chdir("build/")
    # process = subprocess.Popen(
    #     ["ctest"],
    #     stdout=subprocess.PIPE,
    #     stderr=subprocess.PIPE,
    # )
    # stdout, stderr = process.communicate()

    # os.chdir("../../../")