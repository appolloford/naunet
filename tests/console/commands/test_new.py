import tomlkit
from pathlib import Path
from cleo import CommandTester


def test_command_new(tmp_path, application):
    command = application.find("new")
    command_tester = CommandTester(command)

    path = Path(tmp_path / "new_project")
    command_tester.execute(path.as_posix())

    assert "new_project" in (command_tester.io.fetch_output())
    assert "new_project" not in command_tester.io.fetch_error()

    with open(path / "naunet_config.toml") as inpf:
        config = tomlkit.loads(inpf.read())

    general = config["general"]
    chemistry = config["chemistry"]
    chem_element = chemistry["element"]
    chem_species = chemistry["species"]
    assert general["name"] == "new_project"
    assert chem_element["elements"] == []
    assert chem_element["pseudo_elements"] == []
    assert chem_species["allowed"] == []
