# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import os
import re
import sys
import urllib.parse

from .command import Command
from ...species import Species
from ...configuration import Configuration


class InitCommand(Command):
    """
    Initialize Naunet Project directory

    init
        {--name= : Name of the project}
        {--description= : Description of the project}
        {--elements= : List of elements}
        {--pseudo-elements=? : List of pseudo elements}
        {--species= : List of species}
        {--extra-species=? : List of extra required species}
        {--network= : Source of chemical network file}
        {--format= : The format of the chemical network}
        {--dust= : Type of dust model}
        {--heating= : List of heating processes}
        {--cooling= : List of cooling processes}
        {--binding= : List of binding energy of ice species}
        {--yield= : List of photodesorption yields of ice species}
        {--shielding= : List of shielding functions}
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

        config_file = Path.cwd() / "naunet_config.toml"

        if config_file.exists():
            overwrite = self.confirm("Project configure file exists. Overwrite?", False)

            if not overwrite:
                sys.exit()

        name = self.option("name")
        if not name:
            name = Path.cwd().name.lower()

            question = self.create_question(
                "Project name [<comment>{}</comment>]: ".format(name), default=name
            )
            name = self.ask(question)

        description = self.option("description")
        if not description:

            question = self.create_question("Project description :", default="")
            description = self.ask(question)

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

        pseudo_element = self.option("pseudo-elements")
        if pseudo_element:
            pseudo_element = (
                [] if pseudo_element == "null" else pseudo_element.split(",")
            )

        else:
            pseudo_element = Species.default_pseudoelements
            question = self.create_question(
                "Pseudo_elements incuded in the network [<comment>{}</comment>]:".format(
                    pseudo_element
                ),
                default=pseudo_element,
            )
            pseudo_element = self.ask(question)

        species = self.option("species")
        if species:
            species = species.split(",")

        else:
            question = self.create_question(
                "Species incuded in the network [<comment>{}</comment>]:".format([]),
                default=[],
            )
            species = self.ask(question)

        extra_species = self.option("extra-species")
        if extra_species:
            extra_species = [] if extra_species == "null" else extra_species.split(",")

        else:
            question = self.create_question(
                "Extra required species (not in the network) [<comment>{}</comment>]:".format(
                    []
                ),
                default=[],
            )
            extra_species = self.ask(question)

        network = self.option("network")

        if network:
            network = network.split(",")

        else:

            network = [""]

            question = self.create_question(
                "Chemical network file [<comment>{}</comment>]:".format(network),
                default=network,
            )
            network = self.ask(question)

        format = self.option("format")

        if format:
            format = format.split(",")

        else:

            format = [""]

            question = self.create_question(
                "Format of the chemical network [<comment>{}</comment>]:".format(
                    format
                ),
                default=format,
            )
            format = self.ask(question)

        dusttype = self.option("dust")
        if not dusttype:
            dusttype = "none"

        heating = self.option("heating")
        heating = heating.split(",") if heating else []

        cooling = self.option("cooling")
        cooling = cooling.split(",") if cooling else []

        binding = self.option("binding")
        if binding:
            binding = binding.split(",")
            binding = {b.split("=")[0]: b.split("=")[1] for b in binding}
            binding = {s: float(sv) for s, sv in binding.items()}

        yields = self.option("yield")
        if yields:
            yields = yields.split(",")
            yields = {y.split("=")[0]: y.split("=")[1] for y in yields}
            yields = {s: float(sv) for s, sv in yields.items()}

        shielding = self.option("shielding")
        if shielding:
            shielding = shielding.split(",")
            shielding = [it.split(":") for it in shielding]
            shielding = {it[0].strip(): it[1].strip() for it in shielding}

        rate_modifier = self.option("rate-modifier")
        if rate_modifier:
            rate_modifier = [rm.strip() for l in rate_modifier for rm in l.split(",")]
            rate_modifier = [it.split(":") for it in rate_modifier]
            rate_modifier = {it[0].strip(): it[1].strip() for it in rate_modifier}

        ode_modifier = self.option("ode-modifier")
        if ode_modifier:
            ode_modifier = [om.strip() for l in ode_modifier for om in l.split(";")]

        solver = self.option("solver")
        if not solver:
            question = self.create_question(
                "Chemical solver [<comment>{}</comment>]:".format(solver),
                default=solver,
            )
            solver = self.ask(question)

        device = self.option("device")
        if not device:
            question = self.create_question(
                "Computational device [<comment>{}</comment>]:".format(device),
                default=device,
            )
            device = self.ask(question)

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
            device = "gpu"

        config = Configuration(
            name,
            description=description,
            element=element,
            pseudo_element=pseudo_element,
            allowed_species=species,
            required_species=extra_species,
            binding_energy=binding,
            photon_yield=yields,
            network=network,
            format=format,
            heating=heating,
            cooling=cooling,
            shielding=shielding,
            dusttype=dusttype,
            rate_modifier=rate_modifier,
            ode_modifier=ode_modifier,
            solver=solver,
            device=device,
            method=method,
        )

        with open(config_file, "w", encoding="utf-8") as f:
            f.write(config.content)

        render = self.option("render")
        if not render:
            render = self.confirm("Render the source codes now?", False)

        if render:
            self.call("render", "--update-species=false")
