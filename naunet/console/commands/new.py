# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from cleo.helpers import argument
from cleo.helpers import option

from .command import Command
from ...species import Species
from ...configuration import BaseConfiguration


class NewCommand(Command):
    name = "new"
    description = "Create a blank naunet project folder"
    arguments = [argument("path", "Path to create the project.")]
    options = [
        option("name", None, "Project name."),
        option("description", None, "Project description."),
    ]

    def __init__(self):
        super(NewCommand, self).__init__()

    def handle(self):
        from pathlib import Path

        path = Path(self.argument("path"))
        if not path.is_absolute():
            # we do not use resolve here due to compatibility issues
            # for path.resolve(strict=False)
            path = Path.cwd().joinpath(path)

        name = self.option("name") or path.name

        if path.exists() and list(path.glob("*")):
            # Directory is not empty. Aborting.
            raise RuntimeError(
                f"Destination <fg=yellow>{path}</> exists and is not empty"
            )

        path.mkdir(parents=True)

        description = self.option("description") or ""

        config = BaseConfiguration(name, description=description)

        with open(path / "naunet_config.toml", "w", encoding="utf-8") as outf:
            outf.write(config.content)

        self.line(
            f"Created project <info>{config._name}</> in"
            f" <fg=blue>{path.as_posix()}</>"
        )
