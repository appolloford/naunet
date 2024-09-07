# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import ast
import os
import re
import sys
import urllib.parse

from cleo.helpers import argument
from cleo.helpers import option

from .command import Command
from ...species import Species
from ...configuration import BaseConfiguration


class InitCommand(Command):
    name = "init"
    description = "Initialize Naunet Project directory"
    options = [
        option("name", None, "Project name.", flag=False),
        option("description", None, "Project description.", flag=False),
        option(
            "loading",
            None,
            "Python module to be loaded.",
            flag=False,
            value_required=False,
        ),
        option("elements", None, "List of elements.", flag=False, value_required=False),
        option(
            "pseudo-elements",
            None,
            "List of pseudo elements.",
            flag=False,
            value_required=False,
        ),
        option(
            "element-replacement",
            None,
            "Table of elements name to be replaced.",
            flag=False,
            value_required=False,
        ),
        option("surface-prefix", None, "Prefix of surface species.", flag=False),
        option("bulk-prefix", None, "Prefix of bulk species.", flag=False),
        option(
            "allowed-species",
            None,
            "List of allowed species.",
            flag=False,
            value_required=False,
        ),
        option(
            "extra-species",
            None,
            "List of extra required species.",
            flag=False,
            value_required=False,
        ),
        option(
            "binding",
            None,
            "List of ice species binding energies.",
            flag=False,
            value_required=False,
        ),
        option(
            "yield",
            None,
            "List of ice species photondesorption yields.",
            flag=False,
            value_required=False,
        ),
        option("grain-symbol", None, "Dust grain symbol.", flag=False),
        option(
            "grain-model",
            None,
            "Dust grain model type.",
            flag=False,
            value_required=False,
        ),
        option(
            "network-files",
            None,
            "Chemical network files.",
            flag=False,
            value_required=False,
        ),
        option(
            "file-formats",
            None,
            "Chemical network files formats.",
            flag=False,
            value_required=False,
        ),
        option(
            "heating",
            None,
            "Needed heating processes.",
            flag=False,
            value_required=False,
        ),
        option(
            "cooling",
            None,
            "Needed cooling processes.",
            flag=False,
            value_required=False,
        ),
        option(
            "shielding",
            None,
            "Needed shielding functions.",
            flag=False,
            value_required=False,
        ),
        option(
            "rate-modifier",
            None,
            "Reaction rates to be modified.",
            flag=False,
            value_required=False,
            multiple=True,
        ),
        option(
            "ode-modifier",
            None,
            "ODEs to be modified.",
            flag=False,
            value_required=False,
            multiple=True,
        ),
        option("solver", None, "ODE solver.", flag=False),
        option("device", None, "Device.", flag=False),
        option(
            "method",
            None,
            "Linear solver used in CVode or algorithm in ODEInt.",
            flag=False,
        ),
        option("render", None, "Render code after initialization."),
        option("render-force", None, "Force render."),
    ]

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

        loading = self.option("loading")
        loading = self.validate(loading, "Python modules to be loaded", "")
        loading = [l.strip() for l in loading.split(",") if l]

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

        replacement = self.option("element-replacement")
        replacement = self.validate(replacement, "Elements to be replaced", "")
        replacement = replacement.split(",")
        replacement = {r.split(":")[0]: r.split(":")[1] for r in replacement if r}
        replacement = {r.strip(): rv.strip() for r, rv in replacement.items()}

        surface_prefix = self.option("surface-prefix")
        surface_prefix = self.validate(surface_prefix, "Prefix of surface species", "#")

        bulk_prefix = self.option("bulk-prefix")
        bulk_prefix = self.validate(bulk_prefix, "Prefix of bulk species", "@")

        allowed_species = self.option("allowed-species")
        allowed_species = self.validate(
            allowed_species, "Species to be included in the network", ""
        )
        allowed_species = [s.strip() for s in allowed_species.split(",") if s]

        extra_species = self.option("extra-species")
        extra_species = self.validate(
            extra_species,
            "Extra required species (not exist in the network)",
            "",
        )
        extra_species = [e.strip() for e in extra_species.split(",") if e]

        binding = self.option("binding")
        binding = self.validate(binding, "Binding energies", "")
        binding = binding.split(",")
        binding = {b.split("=")[0]: float(b.split("=")[1]) for b in binding if b}

        yields = self.option("yield")
        yields = self.validate(yields, "Photon desorption yields", "")
        yields = yields.split(",")
        yields = {y.split("=")[0]: float(y.split("=")[1]) for y in yields if y}

        grain_symbol = self.option("grain-symbol")
        grain_symbol = self.validate(grain_symbol, "Grain symbol", "GRAIN")

        species_kwargs = {
            "grain_symbol": grain_symbol,
            "surface_prefix": surface_prefix,
            "bulk_prefix": bulk_prefix,
        }

        grain_model = self.option("grain-model")
        grain_model = self.validate(grain_model, "Grain model", "")

        network = self.option("network-files")
        network = self.validate(network, "Chemical network files", "")
        network = [n.strip() for n in network.split(",") if n]

        formats = self.option("file-formats")
        formats = self.validate(formats, "Format of the chemical networks", "")
        formats = [f.strip() for f in formats.split(",") if f]

        heating = self.option("heating")
        heating = self.validate(heating, "Heating models", "")
        heating = [h.strip() for h in heating.split(",") if h]

        cooling = self.option("cooling")
        cooling = self.validate(cooling, "Cooling models", "")
        cooling = [c.strip() for c in cooling.split(",") if c]

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

        ode_modifier_str = self.option("ode-modifier")
        ode_modifier = {}
        for l in ode_modifier_str:
            for om in l.split(";"):
                if not om:
                    break
                key, value = om.split(":")
                fact, rdep = value.split(",")
                rdep = rdep.replace("[", "").replace("]", "").strip().split()
                if ode_modifier.get(key):
                    ode_modifier[key]["factors"].append(fact)
                    ode_modifier[key]["reactants"].append(rdep)
                else:
                    ode_modifier[key] = {
                        "factors": [fact],
                        "reactants": [rdep],
                    }

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

        config = BaseConfiguration(
            name,
            description=description,
            load=loading,
            element=element,
            pseudo_element=pseudo_element,
            replacement=replacement,
            allowed_species=allowed_species,
            required_species=extra_species,
            species_kwargs=species_kwargs,
            binding_energy=binding,
            photon_yield=yields,
            filenames=network,
            formats=formats,
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
