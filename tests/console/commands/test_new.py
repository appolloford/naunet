import shutil
import tomlkit
from pathlib import Path
from cleo import Application
from cleo import CommandTester
from naunet.species import Species
from naunet.console.commands.new import NewCommand


def test_command_new(tmp_path):
    application = Application()
    application.add(NewCommand())

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
        chem_species = chemistry["species"]
        assert general["name"] == "new_project"
        assert chem_species["elements"] == Species.default_elements
        assert chem_species["pseudo_elements"] == Species.default_pseudoelements
        assert chem_species["allowed_species"] == []
