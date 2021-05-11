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
        ism_path = os.path.join(example_path, "ism")

        networklist = [
            "deuterium/1",
            "deuterium/2",
            "deuterium/3",
            "deuterium/4",
            "ism/1",
            "ism/2",
            "ism/3",
        ]

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
        pseudo_element = []
        species = []
        binding_energy = {}
        photon_yield = {}
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

        elif "ism" in case:

            from naunet.examples.ism import (
                element as ism_element,
                pseudo_elements as ism_pelement,
                species as ism_species,
                binding_energy as ism_be,
                photon_yield as ism_yield,
            )

            element = ism_element
            pseudo_element = ism_pelement
            species = ism_species
            binding_energy = ism_be
            photon_yield = ism_yield
            network_source = ism_path

        # copy network file and test
        for file in os.listdir(network_source):

            # skip __init__.py and __pycache__
            if not file.startswith("__"):

                src = os.path.join(network_source, file)
                dest = os.path.join(destination_path, file)

                if os.path.isdir(src):
                    shutil.copytree(src, dest, dirs_exist_ok=True)

                elif os.path.isfile(src):
                    shutil.copyfile(src, dest)

        name = destination if destination else Path.cwd().name.lower()

        # TODO: improve the way to create project in a new directory
        os.chdir(destination_path)

        if case == networklist[0]:

            option = f"--name={name} --description=example"
            option += f" --elements={','.join(element)} --pseudo-elements=o,p,m"
            option += f" --species={','.join(species)}"
            option += " --network=deuterium.krome --database=krome"
            option += " --solver=cvode --device=cpu --method=dense --render"

            self.call("init", option)

        elif case == networklist[1]:

            option = f"--name={name} --description=example"
            option += f" --elements={','.join(element)} --pseudo-elements=o,p,m"
            option += f" --species={','.join(species)}"
            option += " --network=deuterium.krome --database=krome"
            option += " --solver=cvode --device=cpu --method=sparse --render"

            self.call("init", option)

        elif case == networklist[2]:

            option = f"--name={name} --description=example"
            option += f" --elements={','.join(element)} --pseudo-elements=o,p,m"
            option += f" --species={','.join(species)}"
            option += " --network=deuterium.krome --database=krome"
            option += " --solver=cvode --device=cpu --method=cusparse --render"

            self.call("init", option)

        elif case == networklist[3]:

            option = f"--name={name} --description=example"
            option += f" --elements={','.join(element)} --pseudo-elements=o,p,m"
            option += f" --species={','.join(species)}"
            option += " --network=deuterium.krome --database=krome"
            option += " --solver=odeint --device=cpu --method=rosenbrock4 --render"

            self.call("init", option)

        elif case == networklist[4]:

            option = f"--name={name} --description=example"
            option += f" --elements={','.join(element)}"
            option += f" --pseudo-elements={','.join(pseudo_element)}"
            option += f" --species={','.join(species)}"
            option += " --network=rate12_complex.rates --database=leeds"
            option += " --dust=uniform"
            option += f" --binding={','.join(f'{s}={sv}' for s, sv in binding_energy.items())}"
            option += (
                f" --yield={','.join(f'{s}={sv}' for s, sv in photon_yield.items())}"
            )
            option += " --rate-modifier='k[8273] = 0.0'"
            option += " --ode-modifier='double garea = (4*pi*rG*rG) * (y[IDX_GRAIN0I]+y[IDX_GRAINM])'"
            option += " --ode-modifier='double stick1 = (1.0 / (1.0 + 4.2e-2*sqrt(Tgas+Tdust) + 2.3e-3*Tgas - 1.3e-7*Tgas*Tgas))'"
            option += " --ode-modifier='double stick2 = exp(-1741.0/Tgas) / (1.0 + 5e-2*sqrt(Tgas+Tdust) + 1e-14*pow(Tgas, 4.0))'"
            option += " --ode-modifier='double stick = stick1 + stick2'"
            option += " --ode-modifier='double hloss = stick * garea/4.0 * sqrt(8.0*kerg*Tgas/(pi*amu))'"
            option += " --ode-modifier='ydot[IDX_H2I] += 0.5*hloss*y[IDX_HI]; ydot[IDX_HI] -= hloss*y[IDX_HI]'"
            option += " --solver=cvode --device=cpu --method=dense --render"

            self.call("init", option)

        elif case == networklist[5]:

            option = f"--name={name} --description=example"
            option += f" --elements={','.join(element)}"
            option += f" --pseudo-elements={','.join(pseudo_element)}"
            option += f" --species={','.join(species)}"
            option += " --network=rate12_complex.rates --database=leeds"
            option += " --dust=uniform"
            option += f" --binding={','.join(f'{s}={sv}' for s, sv in binding_energy.items())}"
            option += (
                f" --yield={','.join(f'{s}={sv}' for s, sv in photon_yield.items())}"
            )
            option += " --rate-modifier='k[8273] = 0.0'"
            option += " --ode-modifier='double garea = (4*pi*rG*rG) * (y[IDX_GRAIN0I]+y[IDX_GRAINM])'"
            option += " --ode-modifier='double stick1 = (1.0 / (1.0 + 4.2e-2*sqrt(Tgas+Tdust) + 2.3e-3*Tgas - 1.3e-7*Tgas*Tgas))'"
            option += " --ode-modifier='double stick2 = exp(-1741.0/Tgas) / (1.0 + 5e-2*sqrt(Tgas+Tdust) + 1e-14*pow(Tgas, 4.0))'"
            option += " --ode-modifier='double stick = stick1 + stick2'"
            option += " --ode-modifier='double hloss = stick * garea/4.0 * sqrt(8.0*kerg*Tgas/(pi*amu))'"
            option += " --ode-modifier='ydot[IDX_H2I] += 0.5*hloss*y[IDX_HI]; ydot[IDX_HI] -= hloss*y[IDX_HI]'"
            option += " --solver=cvode --device=cpu --method=sparse --render"

            self.call("init", option)

        elif case == networklist[6]:

            option = f"--name={name} --description=example"
            option += f" --elements={','.join(element)}"
            option += f" --pseudo-elements={','.join(pseudo_element)}"
            option += f" --species={','.join(species)}"
            option += " --network=rate12_complex.rates --database=leeds"
            option += " --dust=uniform"
            option += f" --binding={','.join(f'{s}={sv}' for s, sv in binding_energy.items())}"
            option += (
                f" --yield={','.join(f'{s}={sv}' for s, sv in photon_yield.items())}"
            )
            option += " --rate-modifier='k[8273] = 0.0'"
            option += " --ode-modifier='double garea = (4*pi*rG*rG) * (y[IDX_GRAIN0I]+y[IDX_GRAINM])'"
            option += " --ode-modifier='double stick1 = (1.0 / (1.0 + 4.2e-2*sqrt(Tgas+Tdust) + 2.3e-3*Tgas - 1.3e-7*Tgas*Tgas))'"
            option += " --ode-modifier='double stick2 = exp(-1741.0/Tgas) / (1.0 + 5e-2*sqrt(Tgas+Tdust) + 1e-14*pow(Tgas, 4.0))'"
            option += " --ode-modifier='double stick = stick1 + stick2'"
            option += " --ode-modifier='double hloss = stick * garea/4.0 * sqrt(8.0*kerg*Tgas/(pi*amu))'"
            option += " --ode-modifier='ydot[IDX_H2I] += 0.5*hloss*y[IDX_HI]; ydot[IDX_HI] -= hloss*y[IDX_HI]'"
            option += " --solver=cvode --device=cpu --method=cusparse --render"

            self.call("init", option)
