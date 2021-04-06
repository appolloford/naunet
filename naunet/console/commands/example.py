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
    """

    def __init__(self):
        super(ExampleCommand, self).__init__()

    def handle(self):

        from pathlib import Path
        import naunet

        src_parent_path = Path(naunet.__file__).parent
        example_path = os.path.join(src_parent_path, "examples")

        networklist = ["deuterium/1", "deuterium/2", "deuterium/3", "deuterium/4"]

        case = self.choice("Choose an example network", networklist, 0)

        name = Path.cwd().name.lower()
        element = ["e", "H", "D", "He", "C", "N", "O", "GRAIN"]
        species = [
            "C",
            "C+",
            "C-",
            "C2",
            "C2+",
            "C2D",
            "C2D+",
            "C2H",
            "C2H+",
            "C2N",
            "C2N+",
            "C2O+",
            "C3",
            "C3+",
            "CCO",
            "CD",
            "CD+",
            "CD2",
            "CD2+",
            "CH",
            "CH+",
            "CH2",
            "CH2+",
            "CHD",
            "CHD+",
            "CN",
            "CN+",
            "CN-",
            "CNC+",
            "CO",
            "CO+",
            "CO2",
            "CO2+",
            "D",
            "D+",
            "D-",
            "D2O",
            "D2O+",
            "D3O+",
            "DCN",
            "DCN+",
            "DCO",
            "DCO+",
            "DNC",
            "DNC+",
            "DNO",
            "DNO+",
            "DOC+",
            "GRAIN-",
            "GRAIN0",
            "H",
            "H+",
            "H-",
            "H2DO+",
            "H2O",
            "H2O+",
            "H3O+",
            "HCN",
            "HCN+",
            "HCO",
            "HCO+",
            "HD",
            "HD+",
            "HD2O+",
            "HDO",
            "HDO+",
            "HNC",
            "HNC+",
            "HNO",
            "HNO+",
            "HOC+",
            "He",
            "He+",
            "HeD+",
            "HeH+",
            "N",
            "N+",
            "N2",
            "N2+",
            "N2D+",
            "N2H+",
            "N2O",
            "NCO+",
            "ND",
            "ND+",
            "ND2",
            "ND2+",
            "NH",
            "NH+",
            "NH2",
            "NH2+",
            "NHD",
            "NHD+",
            "NO",
            "NO+",
            "NO2",
            "NO2+",
            "O",
            "O+",
            "O-",
            "O2",
            "O2+",
            "O2D",
            "O2D+",
            "O2H",
            "O2H+",
            "OCN",
            "OD",
            "OD+",
            "OD-",
            "OH",
            "OH+",
            "OH-",
            "e-",
            "mD3+",
            "oD2",
            "oD2+",
            "oD2H+",
            "oD3+",
            "oH2",
            "oH2+",
            "oH2D+",
            "oH3+",
            "pD2",
            "pD2+",
            "pD2H+",
            "pD3+",
            "pH2",
            "pH2+",
            "pH2D+",
            "pH3+",
        ]
        if case == networklist[0]:

            netfile = os.path.join(example_path, "deuterium.krome")
            dest = os.path.join(Path.cwd(), "deuterium.krome")
            shutil.copyfile(netfile, dest)

            self.call(
                "init",
                f"--name={name} --description=example --elements={','.join(element)} --pseudo-elements=o,p,m --species={','.join(species)} --network=deuterium.krome --database=krome --solver=cvode --device=cpu --method=dense",
            )

        elif case == networklist[1]:

            netfile = os.path.join(example_path, "deuterium.krome")
            dest = os.path.join(Path.cwd(), "deuterium.krome")
            shutil.copyfile(netfile, dest)

            self.call(
                "init",
                f"--name={name} --description=example --elements={','.join(element)} --pseudo-elements=o,p,m --species={','.join(species)} --network=deuterium.krome --database=krome --solver=cvode --device=cpu --method=sparse",
            )

        elif case == networklist[2]:

            netfile = os.path.join(example_path, "deuterium.krome")
            dest = os.path.join(Path.cwd(), "deuterium.krome")
            shutil.copyfile(netfile, dest)

            self.call(
                "init",
                f"--name={name} --description=example --elements={','.join(element)} --pseudo-elements=o,p,m --species={','.join(species)} --network=deuterium.krome --database=krome --solver=cvode --device=gpu --method=cusparse",
            )

        elif case == networklist[3]:

            netfile = os.path.join(example_path, "deuterium.krome")
            dest = os.path.join(Path.cwd(), "deuterium.krome")
            shutil.copyfile(netfile, dest)

            self.call(
                "init",
                f"--name={name} --description=example --elements={','.join(element)} --pseudo-elements=o,p,m --species={','.join(species)} --network=deuterium.krome --database=krome --solver=odeint --device=cpu --method=rosenbrock4",
            )
