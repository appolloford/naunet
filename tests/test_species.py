import pytest
from naunet.species import Species


def test_known_elements_pseudoelements():
    Species.reset()

    assert Species.known_elements() == Species.default_elements
    assert Species.known_pseudoelements() == Species.default_pseudoelements

    Species("H")

    # change the element list
    elements = ["H", "C", "O"]
    pelements = ["o", "p", "CR"]
    Species.set_known_elements(elements)
    Species.set_known_pseudoelements(pelements)
    assert Species.known_elements() == elements
    assert Species.known_pseudoelements() == pelements

    # add elements
    Species.add_known_elements(["N"])
    assert Species.known_elements() == elements + ["N"]
    Species.add_known_pseudoelements(["PHOTON"])
    assert Species.known_pseudoelements() == pelements + ["PHOTON"]

    # remove elements
    Species.remove_known_elements(["H"])
    assert Species.known_elements() == ["C", "O", "N"]
    Species.remove_known_pseudoelements(["CR"])
    assert Species.known_pseudoelements() == ["o", "p", "PHOTON"]
    Species.reset()


def test_replace_name():
    Species.set_known_elements(["H", "HE", "E"])
    Species._replacement = {"E": "e-", "HE": "He"}
    assert Species("HE").name == "He"
    assert Species("HEH").name == "HeH"
    assert Species("E").name == "e-"
    assert Species("HE").A == 4.0


@pytest.mark.parametrize(
    "name, atom, surface, alias, basename, charge, custalias, eb, mass, mnum, natoms",
    [
        ("e-", False, False, "eM", "e", -1, "Electron", None, 0.0, 0.0, 0),
        ("E", False, False, "EM", "E", -1, "Electron", None, 0.0, 0.0, 0),
        ("H", True, False, "HI", "H", 0, "H", None, 1.007, 1.0, 1),
        ("H+", False, False, "HII", "H", 1, "Hj", None, 1.007, 1.0, 1),
        ("H-", False, False, "HM", "H", -1, "Hk", None, 1.007, 1.0, 1),
        ("H2", False, False, "H2I", "H2", 0, "HH", None, 2.014, 2, 2),
        ("D", True, False, "DI", "D", 0, "D", None, 2.014, 2.0, 1),
        ("oH2D+", False, False, "oH2DII", "oH2D", 1, "H2Dj_ortho", None, 4.028, 4.0, 3),
        ("He", True, False, "HeI", "He", 0, "He", None, 4.002, 4.0, 1),
        ("C", True, False, "CI", "C", 0, "C", None, 12.011, 12.0, 1),
        ("CO", False, False, "COI", "CO", 0, "CO", None, 28.01, 28.0, 2),
        ("N2D+", False, False, "N2DII", "N2D", 1, "N2Dj_ortho", None, 30.028, 30, 3),
        ("Si++++", False, False, "SiIIIII", "Si", 4, "Sijjjj", None, 28.086, 28.0, 1),
        ("#H", False, True, "GHI", "H", 0, "GH", 600.0, 1.007, 1.0, 1),
        ("#CO", False, True, "GCOI", "CO", 0, "GCO", 1150.0, 28.01, 28.0, 2),
        ("#CH4", False, True, "GCH4I", "CH4", 0, "GCH4", 1090.0, 16.039, 16.0, 5),
        # Grain is atom in case renorm number density is required
        # It has no default mass or massnumber, usually they are not needed
        ("GRAIN0", True, False, "GRAIN0I", "GRAIN0", 0, "grain", None, 0.0, 0.0, 1),
        ("GRAIN-", False, False, "GRAINM", "GRAIN", -1, "grain-", None, 0.0, 0.0, 1),
    ],
)
def test_properties_default(
    name, atom, surface, alias, basename, charge, custalias, eb, mass, mnum, natoms
):
    Species.reset()
    species = Species(name)
    assert species.alias == alias
    species.alias = custalias
    assert species.alias == custalias
    assert species.basename == basename
    assert species.charge == charge
    if eb is None:
        with pytest.raises(ValueError, match="not ice species"):
            species.binding_energy
    else:
        assert species.binding_energy == eb
        assert species.eb == eb  # test alias
    assert f"{species.mass:.3f}" == f"{mass:.3f}"
    assert species.massnumber == mnum
    assert species.A == mnum  # test alias
    assert species.n_atoms == natoms
    assert species.is_atom == atom
    assert species.is_surface == surface
    assert f"{species:<12}" == f"{species.name:<12}"


@pytest.mark.parametrize(
    "name, atom, surface, alias, basename, charge, custalias, eb, mass, mnum",
    [
        ("GH", False, True, "GHI", "H", 0, "GH", 600.0, 1.007, 1.0),
        ("GCO", False, True, "GCOI", "CO", 0, "GCO", 1150.0, 28.01, 28.0),
        ("GCH4", False, True, "GCH4I", "CH4", 0, "GCH4", 1090.0, 16.039, 16.0),
    ],
)
def test_properties_surface_prefix_G(
    name, atom, surface, alias, basename, charge, custalias, eb, mass, mnum
):
    species = Species(name, surface_prefix="G")
    assert species.alias == alias
    species.alias = custalias
    assert species.alias == custalias
    assert species.basename == basename
    assert species.charge == charge
    if eb is None:
        with pytest.raises(ValueError, match="not ice species"):
            species.binding_energy
    else:
        assert species.binding_energy == eb
        assert species.eb == eb  # test alias
    assert f"{species.mass:.3f}" == f"{mass:.3f}"
    assert species.massnumber == mnum
    assert species.A == mnum  # test alias
    assert species.is_atom == atom
    assert species.is_surface == surface
    assert f"{species:<12}" == f"{species.name:<12}"


def test_element_counting():
    Species.reset()
    assert Species("H2").element_count.get("H") == 2
    assert Species("N2D+").element_count.get("N") == 2
    assert Species("N2D+").element_count.get("D") == 1
    assert Species("C10").element_count.get("C") == 10
    assert Species("CO2").element_count.get("O") == 2
    Species.add_known_elements(["13C"])
    assert Species("H213CO").element_count.get("13C") == 1
    assert Species("H213CO").element_count.get("H") == 2
    # special case
    assert Species("GRAIN0").element_count.get("GRAIN") == 1
    assert Species("GRAIN-").element_count.get("GRAIN") == 1


def test_enthalpy():
    Species.reset()
    assert Species("#H2O").enthalpy == -238.9
    assert Species("#CH").enthalpy == 592.5


@pytest.mark.parametrize(
    "name1, gsym1, surf1, buik1, name2, gsym2, surf2, buik2",
    [
        ("E", "GRAIN", "#", "@", "e-", "GRAIN", "#", "@"),
        ("#C", "GRAIN", "#", "@", "GC", "GRAIN", "G", "@"),
    ],
)
def test_eq(name1, gsym1, surf1, buik1, name2, gsym2, surf2, buik2):
    spec1 = Species(name1, gsym1, surf1, buik1)
    spec2 = Species(name2, gsym2, surf2, buik2)
    assert spec1 == spec2
    assert spec1 in [spec2]
    assert set([spec1]) == set([spec2])
    assert spec1.is_grain == spec2.is_grain
    assert spec1.is_surface == spec2.is_surface
    Species.add_known_elements(["Fake"])
    assert [Species("Fake"), spec1].index(spec2) == 1
    Species.reset()
