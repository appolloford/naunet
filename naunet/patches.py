from __future__ import annotations
import logging
from copy import copy
from dataclasses import dataclass
from jinja2 import Environment, PackageLoader, Template
from pathlib import Path
from typing import TYPE_CHECKING

from .species import Species
from .reactions import Reaction
from .reactiontype import ReactionType
from .templateloader import NetworkInfo
from .utilities import _prefix, _suffix, _stmwrap

if TYPE_CHECKING:
    from .network import Network


def patch_factory(target: str, device: str, source: str = None):

    if target.lower() == "enzo":
        return EnzoPatch(device, source)
    else:
        raise ValueError("Not supported target!")


class Patch:
    def __init__(self, device: str, source: str) -> None:

        loader = PackageLoader("naunet", source)
        self._env = Environment(loader=loader)
        self._env.filters["stmwrap"] = _stmwrap
        self._env.filters["prefix"] = _prefix
        self._env.filters["suffix"] = _suffix
        self._env.globals.update(zip=zip)
        self._env.trim_blocks = True
        self._env.rstrip_blocks = True

        self._device = device

    @property
    def templates(self):
        templates = self._env.list_templates()
        templates.append("derived_fields.py")
        return templates

    def _render_derived_field(self, info: NetworkInfo, path: Path | str = "./") -> None:
        """
        Generate the derived fields of number density of species for yt. Saved
        in name of "derived_fields_of_network.py"

        Args:
            path (Path | str, optional): Path to save output file. Defaults to "./".

        """

        derived_field_map = []
        derived_species_field = []
        species = info.species

        for s in species:
            alias = "Electron" if s.is_electron else s.alias
            derived_field_map.append(f"    '{s.name}': '{alias}_ndensity',")
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
                        f"    num_unit = dunit / mh_cgs / {1.0 if s.is_electron else s.A}",
                        f"    arr = (num_unit*data['{alias}_Density']).to_ndarray()",
                        f"    return arr",
                    ]
                )
            )

        for ele in info.elements:
            specalias = ["Electron" if s.is_electron else s.alias for s in species]
            specnatom = [s.element_count.get(ele.name, 0) for s in species]
            eleabund = [
                f"{natom}*data['{alias}_ndensity']"
                for natom, alias in zip(specnatom, specalias)
                if natom
            ]
            eleabundstr = " + ".join(["0.0", *eleabund])
            eleabundstr = _stmwrap(eleabundstr, 70, 10)
            derived_field_map.append(f"    'Elem{ele}': 'element_{ele}_ndensity',")
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

            # element abundance on surface
            speconsurface = [s.is_surface for s in species]
            iceeleabund = [
                f"{natom}*data['{alias}_ndensity']"
                for natom, surf, alias in zip(specnatom, speconsurface, specalias)
                if natom and surf
            ]
            iceeleabundstr = " + ".join(["0.0", *iceeleabund])
            iceeleabundstr = _stmwrap(iceeleabundstr, 70, 10)
            derived_field_map.append(
                f"    'IceElem{ele}': 'surface_element_{ele}_ndensity',"
            )
            derived_species_field.append(
                "\n".join(
                    [
                        f"@derived_field(name='surface_element_{ele}_ndensity', sampling_type='cell')",
                        f"def surface_element_{ele}_ndensity(field, data):",
                        f"    if 'enzo' not in data.ds.dataset_type:",
                        f"        return",
                        f"    if data.ds.parameters['MultiSpecies'] < 4:",
                        f"        return",
                        f"    arr = ({iceeleabundstr})",
                        f"    return arr",
                    ]
                )
            )

        derived_field_map_str = "\n".join(
            ["derived_field_map = {", *derived_field_map, "}"]
        )
        with open(path / "derived_fields_of_network.py", "w") as outf:
            outf.write("import yt\n")
            outf.write("from yt import derived_field\n\n")
            outf.write("mh_cgs = float(yt.units.mh_cgs)\n\n")
            outf.write(derived_field_map_str)
            outf.write("\n\n")
            outf.write("\n\n".join(derived_species_field))


class EnzoPatch(Patch):

    _enzo_required_elements = ["e", "H", "D", "He", "C", "O", "Si"]

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

    def __init__(self, device: str, source: str = None) -> None:
        super().__init__(device, source or "templates/patches/enzo")

    def _render(
        self,
        template: Template,
        info: NetworkInfo = None,
        species_group: SpeciesGroups = None,
        save: bool = True,
        path: Path | str = None,
    ) -> None:

        result = template.render(
            network=info,
            species=species_group,
            device=self._device,
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
        network: Network,
        templates: list[str] = None,
        save: bool = True,
        path: Path | str = None,
    ) -> None:

        info = NetworkInfo(
            network.elements,
            network.species,
            network.reactions or [Reaction(reaction_type=ReactionType.DUMMY)],
            network.heating,
            network.cooling,
            network.grains,
            network.shielding,
        )

        known_elements = Species.known_elements().copy()
        Species.add_known_elements(list(Species._replacement.values()))
        for e in self._enzo_required_elements:
            if e not in known_elements:
                Species.add_known_elements([e])
                logging.info(
                    f"Temporarily add Enzo required element {e} "
                    "into known element list"
                )

        # duplicate the species in network and make them use the alias in enzo
        species_network = [copy(s) for s in network.species]

        species_enzo = [Species(s) for s in EnzoPatch.enzo_defined_species_name]
        species_grackle = [Species(s) for s in EnzoPatch.grackle_species_name]

        for s, alias in zip(species_grackle, self.grackle_defined_alias):
            s.alias = alias
            for s1 in species_network:
                if s == s1:
                    s1.alias = alias

        species_intersect_enzo = [s for s in species_network if s in species_enzo]
        species_intersect_grackle = [s for s in species_network if s in species_grackle]
        species_diff_enzo = [s for s in species_network if s not in species_enzo]
        species_diff_grackle = [s for s in species_network if s not in species_grackle]

        species_group = self.SpeciesGroups(
            species_enzo,
            species_grackle,
            species_network,
            species_intersect_enzo,
            species_intersect_grackle,
            species_diff_enzo,
            species_diff_grackle,
        )

        Species.set_known_elements(known_elements)

        templates = templates or self.templates

        for tmplname in templates:
            if tmplname == "derived_fields.py":
                self._render_derived_field(info, path)
            else:
                tmpl = self._env.get_template(tmplname)
                self._render(tmpl, info, species_group, save, path)
