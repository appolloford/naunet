from __future__ import annotations
from copy import copy
from dataclasses import dataclass
from pathlib import Path
from jinja2 import Environment, PackageLoader, Template

from .species import Species
from .templateloader import NetworkInfo
from .utilities import _prefix, _suffix, _stmwrap


class PatchFactory:
    def __init__(
        self,
        netinfo: NetworkInfo,
        target: str,
        device: str,
        cap_species_name: bool = False,
        source: str = None,
    ) -> None:

        if target.lower() == "enzo":
            self.patch = EnzoPatch(netinfo, device, cap_species_name, source)
        else:
            raise ValueError("Not supported target!")

    def render(
        self,
        templates: list[str] = None,
        save: bool = True,
        path: Path | str = "./",
    ):
        self.patch.render(templates, save, path)

    @property
    def templates(self):
        return self.patch.templates


class EnzoPatch:

    # define both e- and E- to exclude electron
    enzo_defined_species_name = [
        "e-",
        "H",
        "H+",
        "He",
        "He+",
        "He++",
        "H-",
        "H2",
        "H2+",
        "D",
        "D+",
        "HD",
        "C",
        "C+",
        "O",
        "O+",
        "Si",
        "Si+",
        "Si++",
        "CH",
        "CH2",
        "CH3+",
        "C2",
        "CO",
        "HCO+",
        "OH",
        "H2O",
        "O2",
    ]

    # repeated define e- and E- to exclude electron
    grackle_species_name = [
        "e-",
        "H",
        "H+",
        "He",
        "He+",
        "He++",
        "H-",
        "H2",
        "H2+",
        "D",
        "D+",
        "HD",
    ]

    grackle_defined_alias = [
        "De",
        "HI",
        "HII",
        "HeI",
        "HeII",
        "HeIII",
        "HM",
        "H2I",
        "H2II",
        "DI",
        "DII",
        "HDI",
    ]

    @dataclass
    class SpeciesGroups:
        enzo: list[Species]
        grackle: list[Species]
        network: list[Species]
        network_int_enzo: list[Species]
        network_int_grackle: list[Species]
        network_diff_enzo: list[Species]
        network_diff_grackle: list[Species]

    def __init__(
        self,
        netinfo: NetworkInfo,
        device: str,
        cap_species_name: bool = False,
        source: str = None,
    ) -> None:

        source = source or "templates/patches/enzo"
        loader = PackageLoader("naunet", source)
        self._env = Environment(loader=loader)
        self._env.filters["stmwrap"] = _stmwrap
        self._env.filters["prefix"] = _prefix
        self._env.filters["suffix"] = _suffix
        self._env.globals.update(zip=zip)
        self._env.trim_blocks = True
        self._env.rstrip_blocks = True

        self._network_info = netinfo
        self.device = device

        enzo_species = [
            Species(s.upper() if cap_species_name else s)
            for s in EnzoPatch.enzo_defined_species_name
        ]
        grackle_species = [
            Species(s.upper() if cap_species_name else s)
            for s in EnzoPatch.grackle_species_name
        ]
        # duplicate the species in network and make them use the alias in enzo
        network_species = [copy(s) for s in netinfo.species]

        for s, alias in zip(grackle_species, self.grackle_defined_alias):
            s.alias = alias
            for s1 in network_species:
                if s == s1:
                    s1.alias = alias

        network_int_enzo_species = [s for s in network_species if s in enzo_species]
        network_int_grackle_species = [
            s for s in network_species if s in grackle_species
        ]

        network_diff_enzo_species = [
            s for s in network_species if s not in enzo_species
        ]
        network_diff_grackle_species = [
            s for s in network_species if s not in grackle_species
        ]

        self._species_group = self.SpeciesGroups(
            enzo_species,
            grackle_species,
            network_species,
            network_int_enzo_species,
            network_int_grackle_species,
            network_diff_enzo_species,
            network_diff_grackle_species,
        )

    def _render(
        self,
        template: Template,
        save: bool = True,
        path: Path | str = None,
    ) -> None:

        result = template.render(
            network=self._network_info,
            species=self._species_group,
            device=self.device,
        )
        name = template.name.replace(".j2", "")

        if save:
            path = Path(path)
            hydro_rk_path = path / "hydro_rk"

            for p in [path, hydro_rk_path]:
                if not p.exists():
                    p.mkdir(parents=True)

            print(path / name)
            with open(path / name, "w") as outf:
                outf.write(result)
        else:
            print(result)

    def render(
        self,
        templates: list[str] = None,
        save: bool = True,
        path: Path | str = None,
    ) -> None:

        templates = templates or self.templates
        print(templates)

        for tmplname in templates:
            if tmplname == "derived_fields.py":
                self._render_derived_field(path)
            else:
                tmpl = self._env.get_template(tmplname)
                self._render(tmpl, save, path)

    @property
    def templates(self):
        templates = self._env.list_templates()
        templates.append("derived_fields.py")
        return templates

    def _render_derived_field(self, path: Path | str = "./") -> None:
        """
        Generate the derived fields of number density of species for yt. Saved
        in name of "derived_fields_of_network.py"

        Args:
            path (Path | str, optional): Path to save output file. Defaults to "./".

        """

        derived_species_field = []
        species = self._network_info.species

        for s in species:
            alias = "Electron" if s.iselectron else s.alias
            derived_species_field.append(
                "\n".join(
                    [
                        f"@derived_field(name='{alias}_ndensity', sampling_type='cell')",
                        f"def {alias}_ndensity(field, data):",
                        f"    if 'enzo' not in data.ds.dataset_type:",
                        f"        return",
                        f"    if data.ds.parameters['MultiSpecies'] < 4:",
                        f"        return",
                        f"    dunit = data.ds.mass_unit/data.ds.length_unit**3",
                        f"    num_unit = dunit / mh_cgs / {1.0 if s.iselectron else s.A}",
                        f"    arr = (num_unit*data['{alias}_Density']).to_ndarray()",
                        f"    return arr",
                    ]
                )
            )

        for ele in Species.known_elements():
            specalias = ["Electron" if s.iselectron else s.alias for s in species]
            specnatom = [s.element_count.get(ele, 0) for s in species]
            eleabund = [
                f"{natom}*data['{alias}_ndensity']"
                for natom, alias in zip(specnatom, specalias)
                if natom
            ]
            eleabundstr = " + ".join(["0.0", *eleabund])
            eleabundstr = _stmwrap(eleabundstr, 70, 10)
            derived_species_field.append(
                "\n".join(
                    [
                        f"@derived_field(name='element_{ele}_ndensity', sampling_type='cell')",
                        f"def element_{ele}_ndensity(field, data):",
                        f"    if 'enzo' not in data.ds.dataset_type:",
                        f"        return",
                        f"    if data.ds.parameters['MultiSpecies'] < 4:",
                        f"        return",
                        f"    arr = ({eleabundstr})",
                        f"    return arr",
                    ]
                )
            )

        with open(path / "derived_fields_of_network.py", "w") as outf:
            outf.write("import yt\n")
            outf.write("from yt import derived_field\n\n")
            outf.write("mh_cgs = float(yt.units.mh_cgs)\n\n")
            outf.write("\n\n".join(derived_species_field))
