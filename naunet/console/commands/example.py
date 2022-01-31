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

from naunet.examples import primordial

from .command import Command


class ExampleCommand(Command):
    """
    Render source codes according to the project setting

    example
        {--show : show the full command only}
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
        fake_path = os.path.join(example_path, "fake")
        primordial_path = os.path.join(example_path, "primordial")
        deuterium_path = os.path.join(example_path, "deuterium")
        ism_path = os.path.join(example_path, "ism")
        cloud_path = os.path.join(example_path, "cloud")

        name = destination if destination else Path.cwd().name.lower()

        networklist = [
            "fake/dense",
            "fake/sparse",
            "fake/cusparse",
            "fake/rosenbrock4",
            "primordial/dense",
            "primordial/sparse",
            "primordial/cusparse",
            "primordial/rosenbrock4",
            "deuterium/dense",
            "deuterium/sparse",
            "deuterium/cusparse",
            "deuterium/rosenbrock4",
            "cloud/dense",
            "cloud/sparse",
            "cloud/rosenbrock4",
            "ism/dense",
            "ism/sparse",
            "ism/cusparse",
        ]

        case_idx = self.option("select")
        if case_idx:
            case = networklist[int(case_idx)]
        else:
            case = self.choice("Choose an example network", networklist, 0)

        element = []
        pseudo_element = []
        species = []
        binding_energy = {}
        photon_yield = {}
        network_source = None

        # Copy the network and the test folder
        if "fake" in case:

            from naunet.examples.fake import (
                element as fake_element,
                species as fake_species,
            )

            element = fake_element
            species = fake_species
            network_source = fake_path
            solver = "odeint" if "rosenbrock4" in case else "cvode"
            device = "gpu" if "cusparse" in case else "cpu"
            method = case.split("/")[-1]

            option = " ".join(
                [
                    f"--name={name} --description=example",
                    f"--elements={','.join(element)} --pseudo-elements=''",
                    f"--species={','.join(species)}",
                    f"--extra-species=''",
                    f"--network=fake.krome --database=krome",
                    f"--solver={solver} --device={device} --method={method}",
                    f"--render",
                ]
            )

        elif "primordial" in case:

            from naunet.examples.primordial import (
                element as primordial_element,
                species as primordial_species,
            )

            element = primordial_element
            species = primordial_species
            cooling = [
                "CIC_HI",
                "CIC_HeI",
                "CIC_HeII",
                "CIC_He_2S",
                "RC_HII",
                "RC_HeI",
                "RC_HeII",
                "RC_HeIII",
                "CEC_HI",
                "CEC_HeI",
                "CEC_HeII",
            ]
            network_source = primordial_path
            solver = "odeint" if "rosenbrock4" in case else "cvode"
            device = "gpu" if "cusparse" in case else "cpu"
            method = case.split("/")[-1]

            option = " ".join(
                [
                    f"--name={name} --description=example",
                    f"--elements={','.join(element)} --pseudo-elements=''",
                    f"--species={','.join(species)}",
                    f"--extra-species=''",
                    f"--cooling={','.join(cooling)}",
                    f"--network=primordial.krome --database=krome",
                    f"--solver={solver} --device={device} --method={method}",
                    f"--render",
                ]
            )

        elif "deuterium" in case:

            from naunet.examples.deuterium import (
                element as deuterium_element,
                species as deuterium_species,
            )

            element = deuterium_element
            species = deuterium_species
            network_source = deuterium_path
            solver = "odeint" if "rosenbrock4" in case else "cvode"
            device = "gpu" if "cusparse" in case else "cpu"
            method = case.split("/")[-1]

            option = " ".join(
                [
                    f"--name={name} --description=example",
                    f"--elements={','.join(element)} --pseudo-elements=o,p,m",
                    f"--species={','.join(species)}",
                    f"--extra-species=''",
                    f"--network=deuterium.krome --database=krome",
                    f"--solver={solver} --device={device} --method={method}",
                    f"--render",
                ]
            )

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
            solver = "odeint" if "rosenbrock4" in case else "cvode"
            device = "gpu" if "cusparse" in case else "cpu"
            method = case.split("/")[-1]

            option = " ".join(
                [
                    f"--name={name} --description=example",
                    f"--elements={','.join(element)}",
                    f"--pseudo-elements={','.join(pseudo_element)}",
                    f"--species={','.join(species)}",
                    f"--extra-species=''",
                    f"--network=rate12_complex.rates --database=leeds",
                    f"--dust=uniform",
                    f"--binding={','.join(f'{s}={sv}' for s, sv in binding_energy.items())}",
                    f"--yield={','.join(f'{s}={sv}' for s, sv in photon_yield.items())}",
                    f"--shielding='H2: L96Table, CO: V09Table, N2: L13Table'",
                    f"--rate-modifier='8274: 0.0'",
                    f"--ode-modifier='double garea = (4*pi*rG*rG) * (y[IDX_GRAIN0I]+y[IDX_GRAINM])'",
                    f"--ode-modifier='double stick1 = (1.0 / (1.0 + 4.2e-2*sqrt(Tgas+Tdust) + 2.3e-3*Tgas - 1.3e-7*Tgas*Tgas))'",
                    f"--ode-modifier='double stick2 = exp(-1741.0/Tgas) / (1.0 + 5e-2*sqrt(Tgas+Tdust) + 1e-14*pow(Tgas, 4.0))'",
                    f"--ode-modifier='double stick = stick1 + stick2'",
                    f"--ode-modifier='double hloss = stick * garea/4.0 * sqrt(8.0*kerg*Tgas/(pi*amu))'",
                    f"--ode-modifier='ydot[IDX_H2I] += 0.5*hloss*y[IDX_HI]; ydot[IDX_HI] -= hloss*y[IDX_HI]'",
                    f"--solver={solver} --device={device} --method={method}",
                    f"--render",
                ]
            )

        elif "cloud" in case:

            from naunet.examples.cloud import (
                elements as cloud_elements,
                pseudo_elements as cloud_pelements,
                species as cloud_species,
                binding_energy as cloud_be,
            )

            elements = cloud_elements
            pseudo_elements = cloud_pelements
            species = cloud_species
            binding_energy = cloud_be
            network_source = cloud_path
            solver = "odeint" if "rosenbrock4" in case else "cvode"
            device = "gpu" if "cusparse" in case else "cpu"
            method = case.split("/")[-1]

            option = " ".join(
                [
                    f"--name={name} --description=example",
                    f"--elements={','.join(elements)}",
                    f"--pseudo-elements={','.join(pseudo_elements)}",
                    f"--species={','.join(species)}",
                    f"--extra-species=''",
                    f"--network=reactions.ucl --database=uclchem",
                    f"--dust=RR07",
                    f"--binding={','.join(f'{s}={sv}' for s, sv in binding_energy.items())}",
                    f"--shielding='CO: VB88Table'",
                    f"--ode-modifier='ydot[IDX_H2I] += H2formation*y[IDX_HI] - H2dissociation*y[IDX_H2I]/nH'",
                    f"--ode-modifier='ydot[IDX_HI] += 2.0*(H2dissociation*y[IDX_H2I]/nH - H2formation*y[IDX_HI])'",
                    f"--solver={solver} --device={device} --method={method}",
                    f"--render",
                ]
            )

        if self.option("show"):
            print("naunet init ", option)

        else:

            # Check whether the test folder exists
            prefix = os.path.join(destination_path, "test")

            if os.path.exists(prefix):
                if not os.path.isdir(prefix):
                    raise FileNotFoundError(
                        errno.ENOENT, os.strerror(errno.ENOENT), prefix
                    )

                elif os.listdir(prefix):
                    overwrite = self.confirm(
                        f"Non-empty test directory. Overwrite?", False
                    )

                    if not overwrite:
                        sys.exit()

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

            # TODO: improve the way to create project in a new directory
            os.chdir(destination_path)

            self.call("init", option)