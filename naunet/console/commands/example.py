# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import errno
import os
import re
import sys
import urllib.parse
import shutil

from .command import Command
from jinja2 import Environment, PackageLoader


class ExampleCommand(Command):
    """
    Render source codes according to the project setting

    example
        {--dry : show the full command only, not create example}
        {--path= : The path to create the project at, create locally if not assigned}
        {--select= : select the example}
        {--render-force : force render}
    """

    def __init__(self):
        super(ExampleCommand, self).__init__()

    def handle(self):

        from pathlib import Path
        import naunet

        path = Path(self.option("path")) if self.option("path") else Path.cwd()
        if not path.is_absolute():
            # we do not use resolve here due to compatibility issues
            # for path.resolve(strict=False)
            path = Path.cwd().joinpath(path)

        if not path.exists():
            path.mkdir(parents=True)

        name = path.name

        networklist = [
            "empty/dense",
            "empty/sparse",
            "empty/cusparse",
            "empty/rosenbrock4",
            "minimal/dense",
            "minimal/sparse",
            "minimal/cusparse",
            "minimal/rosenbrock4",
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

        import importlib

        example = case.split("/")[0]
        examplesrc = naunet_path / "examples" / example
        examplemod = importlib.import_module(f"naunet.examples.{example}")

        description = examplemod.description
        elements = examplemod.elements
        pseudo_elements = examplemod.pseudo_elements
        species = examplemod.species
        extra_species = examplemod.extra_species
        network = examplemod.network
        format = examplemod.format
        heating = examplemod.heating
        cooling = examplemod.cooling
        shielding = examplemod.shielding
        grain_model = examplemod.grain_model
        binding = examplemod.binding_energy
        photon_yield = examplemod.photon_yield
        rate_modifier = examplemod.rate_modifier
        ode_modifier = examplemod.ode_modifier
        solver = "odeint" if "rosenbrock4" in case else "cvode"
        device = "gpu" if "cusparse" in case else "cpu"
        method = case.split("/")[-1]

        shieldingstr = ",".join(f"{key}: {val}" for key, val in shielding.items())
        bindingstr = ",".join(f"{s}={sv}" for s, sv in binding.items())
        yieldstr = ",".join(f"{s}={sv}" for s, sv in photon_yield.items())
        ratemodifierstr = ",".join(f"{r}:{rv}" for r, rv in rate_modifier.items())
        odemodifierstr = ";".join(ode_modifier)

        overwrite = self.option("render-force")

        options = " ".join(
            [
                f"--name={name}",
                f"--description='{description}'",
                f"--elements='{','.join(elements)}'",
                f"--pseudo-elements='{','.join(pseudo_elements)}'",
                f"--species='{','.join(species)}'",
                f"--extra-species='{','.join(extra_species)}'",
                f"--network='{network}'",
                f"--format='{format}'",
                f"--heating='{','.join(heating)}'",
                f"--cooling='{','.join(cooling)}'",
                f"--shielding='{shieldingstr}'",
                f"--grain-model='{grain_model}'",
                f"--binding='{bindingstr}'",
                f"--yield='{yieldstr}'",
                f"--rate-modifier='{ratemodifierstr}'" if ratemodifierstr else "",
                f"--ode-modifier='{odemodifierstr}'" if odemodifierstr else "",
                f"--solver={solver}",
                f"--device={device}",
                f"--method={method}",
                f"--render",
                f"--render-force" if overwrite else "",
            ]
        )

        if self.option("dry"):
            print(f"naunet init {options}")
            return

        # copy network file
        src = examplesrc / network
        dest = path / network
        if network and example != "ism":
            shutil.copyfile(src, dest)

        # create project
        os.chdir(path)
        self.call("init", options)

        # Check whether the test folder exists
        prefix = path / "tests"

        if prefix.exists():
            if not os.path.isdir(prefix):
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), prefix)

            elif os.listdir(prefix):
                overwrite = overwrite or self.confirm(
                    f"Non-empty test directory. Overwrite?", False
                )

                if not overwrite:
                    sys.exit()

        else:
            os.mkdir(prefix)

        # get list of test templates
        testpkgpath = f"templates/tests/{example}"
        testenv = Environment(loader=PackageLoader("naunet", testpkgpath))
        tmplnamelist = testenv.list_templates()

        # render tests
        env = Environment(loader=PackageLoader("naunet"))
        tmplpath = f"tests/{example}"

        for tmplname in tmplnamelist:
            tmpl = env.get_template(f"{tmplpath}/{tmplname}")

            dest = path / "tests" / tmplname.replace(".j2", "")
            with open(dest, "w") as outf:
                outf.write(tmpl.render())
