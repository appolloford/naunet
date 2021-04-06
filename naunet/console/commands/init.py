# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import os
import re
import sys
import urllib.parse

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
        {--elements= : List of elements}
        {--pseudo-elements= : List of pseudo elements}
        {--species= : List of species}
        {--network= : Source of chemical network file}
        {--database= : The database/format of the chemical network}
        {--solver=cvode : ODE solver}
        {--device=cpu : Device}
        {--method=dense : Linear solver used in CVode or algorithm in ODEInt}
    """

    def __init__(self):
        super(InitCommand, self).__init__()

    def handle(self):

        from pathlib import Path

        config_file = os.path.join(Path.cwd(), "naunet_config.toml")

        if os.path.exists(config_file):
            overwrite = self.confirm("Project configure file exists. Overwrite?", False)

            if not overwrite:
                sys.exit()

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

        content.add("general", general)

        chemistry = table()

        element = self.option("elements").split(",")
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

        pseudo_element = self.option("pseudo-elements").split(",")
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

        species = self.option("species").split(",")
        if not species:
            question = self.create_question(
                "Species incuded in the network [<comment>{}</comment>]:".format(
                    species
                ),
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

        odesolver = table()

        solver = self.option("solver")
        if not solver:
            question = self.create_question(
                "Chemical solver [<comment>{}</comment>]:".format(solver),
                default=solver,
            )
            solver = self.ask(question)

        odesolver.add("solver", solver)

        device = self.option("device")
        if not device:
            question = self.create_question(
                "Computational device [<comment>{}</comment>]:".format(device),
                default=device,
            )
            device = self.ask(question)

        odesolver.add("device", device)

        method = self.option("method")
        if not method:
            if solver == "cvode":
                question = self.create_question(
                    "Linear solver used in CVode [<comment>{}</comment>]:".format(
                        method
                    ),
                    default=method,
                )
                method = self.ask(question)
            elif solver == "odeint":
                question = self.create_question(
                    "Algorithm used in ODEInt [<comment>{}</comment>]:".format(method),
                    default=method,
                )
                method = self.ask(question)

        if solver == "cvode" and method == "cusparse":
            odesolver["device"] = "gpu"

        odesolver.add("method", method)

        # required = self.option("required")
        # if not required:
        #     required = self.choice(
        #         "List of required functions.", ["fex", "jac", "jtv"], "0", multiple=True
        #     )

        # odesolver.add("required", required)

        content.add("ODEsolver", odesolver)

        with open(config_file, "w", encoding="utf-8") as f:
            f.write(dumps(content))

        render = self.confirm("Render the source codes now?", False)

        if render:

            self.call("render")
