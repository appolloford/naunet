# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import os
import re
import sys
import urllib.parse
import shutil

import tomlkit

from cleo import option
from datetime import datetime
from tomlkit import nl, comment, table, dumps
from tomlkit.toml_file import TOMLFile

from .command import Command
from ...settings import default_element_list, default_pseudo_element_list


class InitCommand(Command):
    """
    Initialize Naunet Project directory

    init
        {--name= : Name of the project}
        {--description= : Description of the project}
        {--solver=cvode : ODE solver}
        {--device=cpu : Device}
        {--elements=* : List of elements}
        {--pseudo-elements=* : List of pseudo elements}
        {--species=* : List of species}
        {--network= : Source of chemical network file}
        {--database= : The database/format of the chemical network}
    """

    def __init__(self):
        super(InitCommand, self).__init__()

    def handle(self):
        from pathlib import Path

        config_file = os.path.join(Path.cwd(), "naunet_config.toml")

        # current date and time
        now = datetime.now()

        date_string = now.strftime("%d/%m/%Y %H:%M:%S")

        content = tomlkit.document()
        content.add(comment("Naunet config document"))
        # content.add(nl())

        general = table()
        general.add("creation_time", date_string)

        name = self.option("name")
        if not name:
            name = Path.cwd().name.lower()

            question = self.create_question(
                "Project name [<comment>{}</comment>]: ".format(name), default=name
            )
            name = self.ask(question)

        general.add("name", name)

        description = self.option("description")
        if not description:

            question = self.create_question("Project description :", default="")
            description = self.ask(question)

        general.add("description", description)

        solver = self.option("solver")
        question = self.create_question(
            "Chemical solver [<comment>{}</comment>]:".format(solver), default=solver
        )
        solver = self.ask(question)

        general.add("solver", solver)

        device = self.option("device")
        question = self.create_question(
            "Computational device [<comment>{}</comment>]:".format(device),
            default=device,
        )
        solver = self.ask(question)

        general.add("device", device)

        content.add("general", general)

        chemistry = table()

        element = self.option("elements")
        if not element:

            element = default_element_list

            question = self.create_question(
                "Elements incuded in the network [<comment>{}</comment>]:".format(
                    element
                ),
                default=element,
            )
            element = self.ask(question)

        chemistry.add("elements", element)

        pseudo_element = self.option("pseudo-elements")
        if not pseudo_element:

            pseudo_element = default_pseudo_element_list

            question = self.create_question(
                "Pseudo_elements incuded in the network [<comment>{}</comment>]:".format(
                    pseudo_element
                ),
                default=pseudo_element,
            )
            pseudo_element = self.ask(question)

        chemistry.add("pseudo_elements", pseudo_element)

        species = self.option("species")
        question = self.create_question(
            "Species incuded in the network [<comment>{}</comment>]:".format(species),
            default=species,
        )
        species = self.ask(question)

        chemistry.add("species", species)

        network = self.option("network")

        if not network:

            network = ""

            question = self.create_question(
                "Chemical network file [<comment>{}</comment>]:".format(network),
                default=network,
            )
            network = self.ask(question)

        chemistry.add("network", network)

        database = self.option("database")

        if not database:

            database = ""

            question = self.create_question(
                "Source database/format of the chemical network [<comment>{}</comment>]:".format(
                    database
                ),
                default=database,
            )
            database = self.ask(question)

        chemistry.add("database", database)

        content.add("chemistry", chemistry)

        with open(config_file, "w", encoding="utf-8") as f:
            f.write(dumps(content))

        render = self.confirm("Render the source codes now?", False)

        if render:

            self.call("render")

        # import naunet

        # src_parent_path = Path(naunet.__file__).parent
        # csrc_path = str(src_parent_path) + "/cxx_src"

        # dest_path = str(Path.cwd())

        # for file in os.listdir(csrc_path):
        #     src = "/".join([csrc_path, file])
        #     dest = "/".join([dest_path, file])
        #     if os.path.isdir(src):
        #         shutil.copytree(src, dest)
        #     elif os.path.isfile(src):
        #         shutil.copyfile(src, dest)
        #     # else:
        #     #     raise TypeError
