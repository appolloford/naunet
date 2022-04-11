# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from .command import Command
from ...species import Species
from ...configuration import Configuration


class NewCommand(Command):
    """
    Create a new Naunet Project directory

    new
        {path : The path to create the project at}
        {--name= : Name of the project}
        {--description= : Description of the project}
        {--elements= : List of elements}
        {--pseudo-elements= : List of pseudo elements}
        {--species= : List of species}
        {--extra-species= : List of extra required species}
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

        element = self.option("elements")
        element = element.split(",") if element else Species.default_elements

        pseudo_element = self.option("pseudo-elements")
        pseudo_element = (
            pseudo_element.split(",")
            if pseudo_element
            else Species.default_pseudoelements
        )

        species = self.option("species")
        species = species.split(",") if species else []

        extra_species = self.option("extra-species")
        extra_species = extra_species.split(",") if extra_species else []

        network = self.option("network")
        network = network.split(",") if network else []

        format = self.option("format")
        format = format.split(",") if format else []

        dusttype = self.option("dust")
        dusttype = dusttype if dusttype else ""

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
        device = self.option("device")
        method = self.option("method")

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

        content = config.content

        with open(config_file, "w", encoding="utf-8") as outf:
            outf.write(content)

        self.line(
            f"Created project <info>{config._name}</> in"
            f" <fg=blue>{path.as_posix()}</>"
        )
