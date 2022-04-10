# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import errno
import os
import re
import sys
import urllib.parse
import shutil

from .command import Command


class ExampleCommand(Command):
    """
    Render source codes according to the project setting

    example
        {--dry : show the full command only, not create example}
        {--path=. : The path to create the project at, create locally if not assigned}
        {--select= : select the example}
    """

    def __init__(self):
        super(ExampleCommand, self).__init__()

    def handle(self):

        from pathlib import Path
        import naunet

        path = Path(self.option("path"))
        if not path.is_absolute():
            # we do not use resolve here due to compatibility issues
            # for path.resolve(strict=False)
            path = Path.cwd().joinpath(path)

        if not path.exists():
            path.mkdir(parents=True)

        name = path.name

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

        naunet_path = Path(naunet.__file__).parent
        network_source = naunet_path / "examples" / case.split("/")[0]

        # Copy the network and the test folder
        if "fake" in case:

            from naunet.examples.fake import (
                element as fake_element,
                species as fake_species,
            )

            element = fake_element
            species = fake_species
            solver = "odeint" if "rosenbrock4" in case else "cvode"
            device = "gpu" if "cusparse" in case else "cpu"
            method = case.split("/")[-1]

            option = " ".join(
                [
                    f"--name={name} --description=example",
                    f"--elements={','.join(element)} --pseudo-elements=''",
                    f"--species={','.join(species)}",
                    f"--extra-species=''",
                    f"--network=fake.krome --format=krome",
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
                    f"--network=primordial.krome --format=krome",
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
            solver = "odeint" if "rosenbrock4" in case else "cvode"
            device = "gpu" if "cusparse" in case else "cpu"
            method = case.split("/")[-1]

            option = " ".join(
                [
                    f"--name={name} --description=example",
                    f"--elements={','.join(element)} --pseudo-elements=o,p,m",
                    f"--species={','.join(species)}",
                    f"--extra-species=''",
                    f"--network=deuterium.krome --format=krome",
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
                    f"--network=rate12_complex.rates --format=leeds",
                    f"--dust=hh93",
                    f"--binding={','.join(f'{s}={sv}' for s, sv in binding_energy.items())}",
                    f"--yield={','.join(f'{s}={sv}' for s, sv in photon_yield.items())}",
                    f"--shielding='H2: L96Table, CO: V09Table, N2: L13Table'",
                    f"--rate-modifier='8274: 0.0'",
                    # f"--ode-modifier='double garea = (4*pi*rG*rG) * (y[IDX_GRAIN0I]+y[IDX_GRAINM])'",
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
                    f"--network=reactions.ucl --format=uclchem",
                    f"--dust=rr07",
                    f"--binding={','.join(f'{s}={sv}' for s, sv in binding_energy.items())}",
                    f"--shielding='CO: VB88Table'",
                    f"--ode-modifier='ydot[IDX_H2I] += H2formation*y[IDX_HI]*nH - H2dissociation*y[IDX_H2I]'",
                    f"--ode-modifier='ydot[IDX_HI] += 2.0*(H2dissociation*y[IDX_H2I] - H2formation*y[IDX_HI]*nH)'",
                    f"--solver={solver} --device={device} --method={method}",
                    f"--render",
                ]
            )

        if self.option("dry"):
            print(f"naunet init {option}")

        else:

            # Check whether the test folder exists
            prefix = path / "test"

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
                    dest = os.path.join(path, file)

                    if os.path.isdir(src):
                        shutil.copytree(src, dest, dirs_exist_ok=True)

                    elif os.path.isfile(src):
                        shutil.copyfile(src, dest)

            # TODO: improve the way to create project in a new directory
            os.chdir(path)

            self.call("init", option)