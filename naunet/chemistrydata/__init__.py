from __future__ import annotations
import csv
from collections import namedtuple
from pathlib import Path


def _read_csv_table(filename: str, name: str) -> list:
    path = Path(__file__).parent
    with open(path / filename, newline="") as inp:
        reader = csv.reader(filter(lambda row: row[0] != "#", inp))
        item = namedtuple(name, next(reader))
        table = list(map(item._make, reader))
    return table


def _read_binding_energy() -> dict[str, int]:
    path = Path(__file__).parent
    with open(path / "rate12_binding_energy.dat", newline="") as inp:
        binding_energy = {}
        for line in inp.readlines():
            if not line.startswith("#"):
                elem, eb, *other = line.split()
                binding_energy.update({elem: float(eb)})
    return binding_energy


# Load the periodic table as a list of namedtuple
# Source: https://gist.github.com/GoodmanSciences/c2dd862cd38f21b0ad36b8f96b4bf1ee
periodic_table = _read_csv_table("periodictable.csv", "Element")
# reference: http://moltensalt.org/references/static/downloads/pdf/stable-isotopes.pdf
isotopes_table = _read_csv_table("isotopestable.csv", "Isotope")
# reference: https://cccbdb.nist.gov/hf0k.asp
enthalpy_table = _read_csv_table("enthalpytable.csv", "Enthalpy")
# reference: http://udfa.ajmarkwick.net/downloads/RATE12_binding_energies.dist.txt
rate12_binding_energy = _read_binding_energy()

user_binding_energy = {}

user_photon_yield = {}


def update_binding_energy(eb: dict):
    user_binding_energy.update(eb)


def update_photon_yield(phyield: dict):
    user_photon_yield.update(phyield)
