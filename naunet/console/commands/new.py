# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from .command import Command
from ...configuration import Configuration


class NewCommand(Command):
    """
    Create a new Naunet Project directory

    new
        {path : The path to create the project at}
        {--name=? : Name of the project}
        {--description=? : Description of the project}
        {--solver=cvode : ODE solver}
        {--device=cpu : Device}
        {--method=dense : Linear solver used in CVode or algorithm in ODEInt}
    """

    def __init__(self):
        super(NewCommand, self).__init__()

    def handle(self):

        from pathlib import Path

        path = Path(self.argument("path"))
        if not path.is_absolute():
            # we do not use resolve here due to compatibility issues
            # for path.resolve(strict=False)
            path = Path.cwd().joinpath(path)

        name = self.option("name")
        if not name:
            name = path.name

        if path.exists() and list(path.glob("*")):
            # Directory is not empty. Aborting.
            raise RuntimeError(
                f"Destination <fg=yellow>{path}</> exists and is not empty"
            )

        path.mkdir(parents=True)

        config_file = path / "naunet_config.toml"

        description = self.option("description")
        if not self.option("description"):
            description = ""
        solver = self.option("solver")
        device = self.option("device")
        method = self.option("method")

        config = Configuration(
            name,
            description=description,
            solver=solver,
            device=device,
            method=method,
        )

        content = config.content

        with open(config_file, "w", encoding="utf-8") as outf:
            outf.write(content)

        self.line(
            f"Created project <info>{config._name}</> in"
            f" <fg=blue>{path.as_posix()}</>"
        )
