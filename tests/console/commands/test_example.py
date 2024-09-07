import os
import subprocess
import pytest
import tomlkit
from pathlib import Path
from cleo.testers.command_tester import CommandTester

GITHUB_ACTIONS = os.getenv("GITHUB_ACTIONS") == "true"


@pytest.mark.skipif(GITHUB_ACTIONS, reason="on github")
def test_command_example(tmp_path, application):
    command = application.find("example")
    command_tester = CommandTester(command)

    os.chdir(tmp_path)
    command_tester.execute("--select=0")

    assert (tmp_path / "naunet_config.toml").exists()
    assert (tmp_path / "include").exists()
    assert (tmp_path / "src").exists()
    assert (tmp_path / "tests").exists()

    with open(tmp_path / "naunet_config.toml") as inpf:
        config = tomlkit.loads(inpf.read())

    general = config["general"]
    # chemistry = config["chemistry"]
    # chem_element = chemistry["element"]
    # chem_species = chemistry["species"]
    assert general["name"] == tmp_path.name
    # assert chem_element["elements"] == []
    # assert chem_element["pseudo_elements"] == []
    # assert chem_species["allowed"] == []

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
