import shutil
import tomlkit
from pathlib import Path
from cleo import Application
from cleo import CommandTester
from naunet.console.commands.new import NewCommand


def test_command_new():
    application = Application()
    application.add(NewCommand())

    command = application.find("new")
    command_tester = CommandTester(command)

    path = Path("test/test_output/new_project")
    command_tester.execute(path.as_posix())

    assert "new_project" in (command_tester.io.fetch_output())
    assert "new_project" not in command_tester.io.fetch_error()

    with open(path / "naunet_config.toml") as inpf:
        config = tomlkit.loads(inpf.read())

        general = config["general"]
        chemistry = config["chemistry"]
        assert general["name"] == "new_project"
        assert chemistry["elements"] == []
        assert chemistry["species"] == []

    shutil.rmtree(path)
