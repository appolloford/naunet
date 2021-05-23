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
from ...species import Species


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
        {--dust= : Type of dust model}
        {--binding= : List of binding energy of ice species}
        {--yield= : List of photodesorption yields of ice species}
        {--rate-modifier=* : List of reaction rates to be changed}
        {--ode-modifier=* : List of ODEs to be changed}
        {--solver=cvode : ODE solver}
        {--device=cpu : Device}
        {--method=dense : Linear solver used in CVode or algorithm in ODEInt}
        {--render : render template immediately}
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

        element = self.option("elements")
        if element:
            element = element.split(",")

        else:
            element = Species.default_elements
            question = self.create_question(
                "Elements incuded in the network [<comment>{}</comment>]:".format(
                    element
                ),
                default=element,
            )
            element = self.ask(question)

        chemistry.add("elements", element)

        pseudo_element = self.option("pseudo-elements")
        if pseudo_element:
            pseudo_element = pseudo_element.split(",")

        else:
            pseudo_element = Species.default_pseudoelements
            question = self.create_question(
                "Pseudo_elements incuded in the network [<comment>{}</comment>]:".format(
                    pseudo_element
                ),
                default=pseudo_element,
            )
            pseudo_element = self.ask(question)

        chemistry.add("pseudo_elements", pseudo_element)

        species = self.option("species")
        if species:
            species = species.split(",")

        else:
            question = self.create_question(
                "Species incuded in the network [<comment>{}</comment>]:".format([]),
                default=[],
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

        dust = table()

        dustype = self.option("dust")
        if not dustype:
            dustype = "none"
        dust.add("type", dustype)

        chemistry.add("dust", dust)

        binding = table()

        binding_list = self.option("binding")
        if binding_list:
            binding_list = binding_list.split(",")
            binding_list = {b.split("=")[0]: b.split("=")[1] for b in binding_list}
            binding_list = {s: float(sv) for s, sv in binding_list.items()}
            binding.update(binding_list)

        chemistry.add("binding_energy", binding)

        yields = table()

        yield_list = self.option("yield")
        if yield_list:
            yield_list = yield_list.split(",")
            yield_list = {y.split("=")[0]: y.split("=")[1] for y in yield_list}
            yield_list = {s: float(sv) for s, sv in yield_list.items()}
            yields.update(yield_list)

        chemistry.add("photon_yield", yields)

        rate_modifier = self.option("rate-modifier")
        if rate_modifier:
            rate_modifier = [rm.strip() for l in rate_modifier for rm in l.split(";")]

        chemistry.add("rate_modifier", rate_modifier)

        ode_modifier = self.option("ode-modifier")
        if ode_modifier:
            ode_modifier = [om.strip() for l in ode_modifier for om in l.split(";")]

        chemistry.add("ode_modifier", ode_modifier)

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

        render = self.option("render")
        if not render:
            render = self.confirm("Render the source codes now?", False)

        if render:

            self.call("render", "--update-species=false")
