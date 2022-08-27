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
        {--elements=? : List of elements}
        {--pseudo-elements=? : List of pseudo elements}
        {--species=? : List of species}
        {--extra-species=? : List of extra required species}
        {--network=? : Source of chemical network file}
        {--format=? : The format of the chemical network}
        {--grain-model=? : Type of dust grain model}
        {--heating=? : List of heating processes}
        {--cooling=? : List of cooling processes}
        {--binding=? : List of binding energy of ice species}
        {--yield=? : List of photodesorption yields of ice species}
        {--shielding=? : List of shielding functions}
        {--rate-modifier=* : List of reaction rates to be changed}
        {--ode-modifier=* : List of ODEs to be changed}
        {--solver= : ODE solver}
        {--device= : Device}
        {--method= : Linear solver used in CVode or algorithm in ODEInt}
        {--render : render template immediately}
        {--render-force : force render}
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
        name = self.validate(name, "Project name", Path.cwd().name.lower())

        description = self.option("description")
        description = self.validate(description, "Project description", "")

        element = self.option("elements")
        element = self.validate(
            element,
            "Elements incuded in the network",
            ", ".join(Species.default_elements),
        )
        element = [e.strip() for e in element.split(",") if e]

        pseudo_element = self.option("pseudo-elements")
        pseudo_element = self.validate(
            pseudo_element,
            "Pseudo_elements incuded in the network",
            ", ".join(Species.default_pseudoelements),
        )
        pseudo_element = [p.strip() for p in pseudo_element.split(",") if p]

        species = self.option("species")
        species = self.validate(species, "Species incuded in the network", "")
        species = [s.strip() for s in species.split(",") if s]

        extra_species = self.option("extra-species")
        extra_species = self.validate(
            extra_species,
            "Extra required species (not exist in the network)",
            "",
        )
        extra_species = [e.strip() for e in extra_species.split(",") if e]

        network = self.option("network")
        network = self.validate(network, "Chemical network files", "")
        network = [n.strip() for n in network.split(",") if n]

        format = self.option("format")
        format = self.validate(format, "Format of the chemical networks", "")
        format = [f.strip() for f in format.split(",") if f]

        grain_model = self.option("grain-model")
        grain_model = self.validate(grain_model, "Grain model", "")

        heating = self.option("heating")
        heating = self.validate(heating, "Heating models", "")
        heating = [h.strip() for h in heating.split(",") if h]

        cooling = self.option("cooling")
        cooling = self.validate(cooling, "Cooling models", "")
        cooling = [c.strip() for c in cooling.split(",") if c]

        binding = self.option("binding")
        binding = self.validate(binding, "Binding energies", "")
        binding = binding.split(",")
        binding = {b.split("=")[0]: float(b.split("=")[1]) for b in binding if b}

        yields = self.option("yield")
        yields = self.validate(yields, "Photon desorption yields", "")
        yields = yields.split(",")
        yields = {y.split("=")[0]: float(y.split("=")[1]) for y in yields if y}

        shielding = self.option("shielding")
        shielding = self.validate(shielding, "Self shielding models", "")
        shielding = shielding.split(",")
        shielding = {s.split(":")[0]: s.split(":")[1] for s in shielding if s}
        shielding = {s.strip(): sv.strip() for s, sv in shielding.items()}

        # To receive multiple values, rate modifier and ode modifier will not ask
        # question even if the value is not provided
        rate_modifier = self.option("rate-modifier")
        rate_modifier = [rm.strip() for l in rate_modifier for rm in l.split(",")]
        rate_modifier = [rm.split(":") for rm in rate_modifier]
        rate_modifier = {rm[0].strip(): rm[1].strip() for rm in rate_modifier}

        ode_modifier = self.option("ode-modifier")
        ode_modifier = [om.strip() for l in ode_modifier for om in l.split(";")]

        solver = self.option("solver")
        solver = self.validate(solver, "Differential equation solver", "cvode")

        device = self.option("device")
        device = self.validate(device, "Computational device", "cpu")

        allowed_method = {
            "cvode": {"cpu": ["dense", "sparse"], "gpu": ["cusparse"]},
            "odeint": {"cpu": ["rosenbrock4"], "gpu": [""]},
        }

        choices = allowed_method.get(solver).get(device)
        method = self.option("method")
        if method is None and choices:
            method = self.choice(f"Method in {solver}", choices, 0)

        elif method not in choices:
            raise ValueError(f"Not support {method} in {solver} on {device}")

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
            grain_model=grain_model,
            rate_modifier=rate_modifier,
            ode_modifier=ode_modifier,
            solver=solver,
            device=device,
            method=method,
        )

        with open(config_file, "w", encoding="utf-8") as f:
            f.write(config.content)

        render = self.option("render") or self.confirm(
            "Render the source codes now?", False
        )

        if render:
            overwrite = self.option("render-force")
            self.call("render", "--force" if overwrite else "")

    def option(self, key=None):

        value = super().option(key)
        if isinstance(value, str):
            value = value.replace("null", "")

        return value

    def validate(self, value, question, default):

        if value is None:

            value = default
            question = self.create_question(
                f"{question} [<comment>{default}</comment>]:",
                default=default,
            )
            value = self.ask(question)

        return value
