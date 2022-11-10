from __future__ import annotations

import tomlkit

from datetime import datetime
from typing import TYPE_CHECKING

from .species import Species

if TYPE_CHECKING:
    from .network import Network

NAUNET_CONFIG_DEFAULT = """\
# Naunet config document

[general]
creation_time = ""
name = ""
description = ""
loads = []

[chemistry]
elements = []
pseudo_elements = []
species = []
extra_species = []
network = []
format = []
heating = []
cooling = []

[chemistry.grain]
model = ""

[chemistry.binding_energy]

[chemistry.photon_yield]

[chemistry.shielding]

[chemistry.rate_modifier]

[chemistry.ode_modifier]

[ODEsolver]
solver = ""
device = ""
method = ""

[summary]
num_of_elements = -1
num_of_species = -1
num_of_grains = -1
num_of_gas_species = -1
num_of_ice_species = -1
num_of_grain_species = -1
num_of_reactions = -1
list_of_elements = []
list_of_species = []
list_of_species_alias = []
list_of_gas_species = []
list_of_ice_species = []
list_of_grain_species = []
"""


class Configuration:
    def __init__(
        self,
        project: str,
        description: str = "",
        load: list[str] = None,
        element: list[str] = None,
        pseudo_element: list[str] = None,
        allowed_species: list[str] = None,
        required_species: list[str] = None,
        binding_energy: dict[str, float] = None,
        photon_yield: dict[str, float] = None,
        network: list[str] = None,
        format: list[str] = None,
        heating: list[str] = None,
        cooling: list[str] = None,
        shielding: dict[str, str] = None,
        grain_model: str = "",
        rate_modifier: dict[int, str] = None,
        ode_modifier: dict[str, dict[str, list[str | list[str]]]] = None,
        solver: str = "cvode",
        device: str = "cpu",
        method: str = "dense",
        instance: Network = None,
    ) -> None:

        self._name = project
        self._description = description
        self._load = load or []
        self._element = element if element else Species.known_elements()
        self._pseudoelement = (
            pseudo_element if pseudo_element else Species.known_pseudoelements()
        )
        self._species = allowed_species if allowed_species else []
        self._extraspecies = required_species if required_species else []
        self._bindingenergy = binding_energy if binding_energy else {}
        self._photonyield = photon_yield if photon_yield else {}
        self._network = network if network else []
        self._format = format if format else []
        self._heating = heating if heating else []
        self._cooling = cooling if cooling else []
        self._shielding = shielding if shielding else {}
        self._grain_model = grain_model
        self._ratemodifier = rate_modifier if rate_modifier else {}
        self._odemodifier = ode_modifier or {}
        self._solver = solver
        self._device = device
        self._method = method
        self._instance = instance

    @property
    def content(self) -> str:
        content = tomlkit.loads(NAUNET_CONFIG_DEFAULT)

        general = content["general"]

        general["creation_time"] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        general["name"] = self._name
        general["description"] = self._description
        general["loads"] = self._load

        chemistry = content["chemistry"]
        chemistry["elements"] = self._element
        chemistry["pseudo_elements"] = self._pseudoelement
        chemistry["species"] = self._species
        chemistry["extra_species"] = self._extraspecies
        chemistry["binding_energy"] = self._bindingenergy
        chemistry["photon_yield"] = self._photonyield
        chemistry["network"] = self._network
        chemistry["format"] = self._format
        chemistry["heating"] = self._heating
        chemistry["cooling"] = self._cooling
        chemistry["shielding"] = self._shielding
        chemistry["grain"]["model"] = self._grain_model
        chemistry["rate_modifier"] = self._ratemodifier
        chemistry["ode_modifier"] = self._odemodifier

        odesolver = content["ODEsolver"]
        odesolver["solver"] = self._solver
        odesolver["device"] = self._device
        odesolver["method"] = self._method

        summary = content["summary"]

        instance = self._instance
        if instance:
            grains = instance.grains
            gas_species = [s.name for s in instance.species if not s.is_surface]
            ice_species = [g.name for g in instance.species if g.is_surface]
            grain_species = (
                [s.name for g in grains for s in g.species] if grains else []
            )
            summary["num_of_elements"] = len(instance.elements)
            summary["num_of_species"] = len(instance.species)
            summary["num_of_grains"] = len(instance.grains)
            summary["num_of_gas_species"] = len(gas_species)
            summary["num_of_ice_species"] = len(ice_species)
            summary["num_of_grain_species"] = len(grain_species)
            summary["num_of_reactions"] = len(instance.reactions)
            summary["list_of_elements"] = [x.name for x in instance.elements]
            summary["list_of_species"] = [x.name for x in instance.species]
            summary["list_of_species_alias"] = [x.alias for x in instance.species]
            summary["list_of_gas_species"] = gas_species
            summary["list_of_ice_species"] = ice_species
            summary["list_of_grain_species"] = grain_species

        return tomlkit.dumps(content)
