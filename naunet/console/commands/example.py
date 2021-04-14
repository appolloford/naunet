# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import errno
import os
import re
import sys
import urllib.parse
import shutil

from cleo import argument, option
from tomlkit import dumps
from tomlkit.toml_file import TOMLFile

from .command import Command


class ExampleCommand(Command):
    """
    Render source codes according to the project setting

    example
        {--select= : select the example}
        {--dest= : the destination directory, created if not existed}
    """

    def __init__(self):
        super(ExampleCommand, self).__init__()

    def handle(self):

        from pathlib import Path
        import naunet

        destination_path = Path.cwd()
        destination = self.option("dest")
        if destination:
            Path(destination).mkdir(parents=True)
            destination_path = destination

        src_parent_path = Path(naunet.__file__).parent
        example_path = os.path.join(src_parent_path, "examples")
        deuterium_path = os.path.join(example_path, "deuterium")

        networklist = ["deuterium/1", "deuterium/2", "deuterium/3", "deuterium/4"]

        case_idx = self.option("select")
        if case_idx:
            case = networklist[int(case_idx)]
        else:
            case = self.choice("Choose an example network", networklist, 0)

        # Check whether the test folder exists
        prefix = os.path.join(destination_path, "test")

        if os.path.exists(prefix):
            if not os.path.isdir(prefix):
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), prefix)

            elif os.listdir(prefix):
                overwrite = self.confirm(f"Non-empty test directory. Overwrite?", False)

                if not overwrite:
                    sys.exit()

        element = []
        species = []
        network_source = None

        # Copy the network and the test folder
        if "deuterium" in case:

            from naunet.examples.deuterium import (
                element as deuterium_element,
                species as deuterium_species,
            )

            element = deuterium_element
            species = deuterium_species
            network_source = deuterium_path

        for file in os.listdir(network_source):

            # skip __init__.py and __pycache__
            if not file.startswith("__"):

                src = os.path.join(network_source, file)
                dest = os.path.join(destination_path, file)

                if os.path.isdir(src):
                    shutil.copytree(src, dest)

                elif os.path.isfile(src):
                    shutil.copyfile(src, dest)

        name = destination if destination else Path.cwd().name.lower()

        # TODO: improve the way to create project in a new directory
        os.chdir(destination_path)

        if case == networklist[0]:

            self.call(
                "init",
                f"--name={name} --description=example --elements={','.join(element)} --pseudo-elements=o,p,m --species={','.join(species)} --network=deuterium.krome --database=krome --solver=cvode --device=cpu --method=dense --render",
            )

        elif case == networklist[1]:

            self.call(
                "init",
                f"--name={name} --description=example --elements={','.join(element)} --pseudo-elements=o,p,m --species={','.join(species)} --network=deuterium.krome --database=krome --solver=cvode --device=cpu --method=sparse --render",
            )

        elif case == networklist[2]:

            self.call(
                "init",
                f"--name={name} --description=example --elements={','.join(element)} --pseudo-elements=o,p,m --species={','.join(species)} --network=deuterium.krome --database=krome --solver=cvode --device=gpu --method=cusparse --render",
            )

        elif case == networklist[3]:

            self.call(
                "init",
                f"--name={name} --description=example --elements={','.join(element)} --pseudo-elements=o,p,m --species={','.join(species)} --network=deuterium.krome --database=krome --solver=odeint --device=cpu --method=rosenbrock4 --render",
            )
