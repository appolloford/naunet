# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import errno
import os
import re
import sys
import urllib.parse
import shutil

import tomlkit

from pathlib import Path
from importlib import util

from cleo.helpers import argument
from cleo.helpers import option
from tomlkit.toml_file import TOMLFile

from naunet.templateloader import TemplateLoader
from naunet.patches import patch_factory

from .command import Command


class RenderCommand(Command):
    name = "render"
    description = "Render source codes according to the project setting"
    options = [
        option("force", "f", "Force render, override existing files."),
        option("patch", None, "Create patch files for target code.", flag=False),
        option("patch-source", None, "Patch source directory.", flag=False),
        option("with-pattern", None, "Render Jacobian pattern."),
    ]

    def __init__(self):
        super(RenderCommand, self).__init__()

    def handle(self):
        config = TOMLFile("naunet_config.toml")

        content = config.read()
        general = content["general"]
        name = general["name"]
        loads = general["loads"]

        if loads:
            for l in loads:
                spec = util.spec_from_file_location(l, Path.cwd() / l)
                module = util.module_from_spec(spec)
                spec.loader.exec_module(module)

        chemistry = content["chemistry"]

        chem_symbol = chemistry["symbol"]
        chem_element = chemistry["element"]
        chem_species = chemistry["species"]
        chem_network = chemistry["network"]
        chem_grain = chemistry["grain"]
        chem_thermal = chemistry["thermal"]

        element = chem_element["elements"]
        pseudo_element = chem_element["pseudo_elements"]
        replacement = chem_element["replacement"]
        allowed_species = chem_species["allowed"]
        extra_species = chem_species["required"]
        binding = chem_species["binding_energy"]
        yields = chem_species["photon_yield"]
        species_kwargs = {
            "grain_symbol": chem_symbol["grain"],
            "surface_prefix": chem_symbol["surface"],
            "bulk_prefix": chem_symbol["bulk"],
        }

        files = chem_network["files"]
        formats = chem_network["formats"]

        heating = chem_thermal["heating"]
        cooling = chem_thermal["cooling"]

        shielding = chemistry["shielding"]
        rate_modifier = chemistry["rate_modifier"]
        ode_modifier = chemistry["ode_modifier"]

        grain_model = chem_grain["model"]

        odesolver = content["ODEsolver"]
        solver = odesolver["solver"]
        method = odesolver["method"]
        device = odesolver["device"]
        # required = odesolver["required"]

        import naunet
        from naunet.species import Species
        from naunet.network import Network, supported_reaction_class
        from naunet.chemistrydata import update_binding_energy, update_photon_yield

        Species._replacement = replacement
        Species.set_known_elements(element)
        Species.set_known_pseudoelements(pseudo_element)
        binding = {
            Species(key, **species_kwargs).name: value for key, value in binding.items()
        }
        yields = {
            Species(key, **species_kwargs).name: value for key, value in yields.items()
        }

        update_binding_energy(binding)
        update_photon_yield(yields)

        rate_modifier = {int(key): value for key, value in rate_modifier.items()}
        net = Network(
            filelist=files,
            fileformats=formats,
            elements=element,
            pseudo_elements=pseudo_element,
            allowed_species=allowed_species,
            required_species=extra_species,
            species_kwargs=species_kwargs,
            grain_model=grain_model,
            heating=heating,
            cooling=cooling,
            shielding=shielding,
            rate_modifier=rate_modifier,
            ode_modifier=ode_modifier,
        )

        dupes, dupidx, first = net.find_duplicate_reaction(mode="short")
        print(f"The following {len(first)} reactions appear multiple times:\n")
        print("".join([fst.react_string for fst in first]))
        print(f"They repeatedly appear in the following {len(dupes)} reactions:\n")
        print("".join([dup.react_string for dup in dupes]))

        patchname = self.option("patch")
        source = self.option("patch-source")
        # include/src/tests are not changed when rendering patches
        if patchname:
            patch = patch_factory(patchname, device, source)
            patch.render(net, path=Path.cwd() / patchname)
            return

        # If not creating patch, check whether include, src, python folders exist
        header_prefix = Path.cwd() / "include"
        source_prefix = Path.cwd() / "src"
        python_prefix = Path.cwd() / "python"
        for prefix in [header_prefix, source_prefix, python_prefix]:
            if prefix.exists():
                if not os.path.isdir(prefix):
                    raise FileNotFoundError(
                        errno.ENOENT, os.strerror(errno.ENOENT), prefix
                    )

                elif os.listdir(prefix):
                    overwrite = self.option("force") or self.confirm(
                        f"Non-empty {prefix.name} directory. Overwrite?", False
                    )

                    if not overwrite:
                        sys.exit()

            else:
                os.mkdir(prefix)

        pattern = self.option("with-pattern")
        tl = TemplateLoader(solver=solver, method=method, device=device)
        tl.render(name, net, path=Path.cwd(), jac_pattern=pattern)

        pkgpath = Path(naunet.__file__).parent

        demo = Path.cwd() / "demo.ipynb"
        if not demo.exists():
            shutil.copyfile(pkgpath / "templates/base/demo.ipynb", demo)

        summary = tomlkit.table()
        all_elements = [e.name for e in net.elements]
        all_species = [x.name for x in net.species]
        all_alias = [x.alias for x in net.species]
        gas_species = [s.name for s in net.species if not s.is_surface]
        ice_species = [g.name for g in net.species if g.is_surface]
        grain_species = [s.name for s in net.species if s.is_grain]
        summary["num_of_elements"] = len(net.elements)
        summary["num_of_species"] = len(net.species)
        summary["num_of_grains"] = len(net.grains)
        summary["num_of_gas_species"] = len(gas_species)
        summary["num_of_ice_species"] = len(ice_species)
        summary["num_of_grain_species"] = len(grain_species)
        summary["num_of_reactions"] = len(net.reactions)
        summary["list_of_elements"] = all_elements
        summary["list_of_species"] = all_species
        summary["list_of_species_alias"] = all_alias
        summary["list_of_gas_species"] = gas_species
        summary["list_of_ice_species"] = ice_species
        summary["list_of_grain_species"] = grain_species

        content["summary"] = summary

        config_file = Path.cwd() / "naunet_config.toml"
        with open(config_file, "w", encoding="utf-8") as f:
            f.write(tomlkit.dumps(content))

        # progress = self.progress_bar()
        # progress.finish()
