from naunet.species import Species


def test_known_elements_pseudoelements():
    Species.reset()

    # list should be empty before any initialization
    assert Species.known_elements() == []
    assert Species.known_pseudoelements() == []

    # Use default element list if class is used before set them
    Species("H")
    assert Species.known_elements() == Species.default_elements
    assert Species.known_pseudoelements() == Species.default_pseudoelements

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


def test_alias():
    Species.reset()
    assert Species("H2").alias == "H2I"
    assert Species("CO").alias == "COI"
    assert Species("e-").alias == "eM"
    assert Species("H-").alias == "HM"
    assert Species("H+").alias == "HII"
    assert Species("oH2D+").alias == "oH2DII"
    assert Species("N2D+").alias == "N2DII"
    assert Species("#CO").alias == "GCOI"


def test_alias_setter():
    Species.reset()
    spec = Species("H2")
    spec.alias = "HH"
    assert spec.alias == "HH"


def test_basename():
    Species.reset()
    assert Species("H2").basename == "H2"
    assert Species("CO").basename == "CO"
    assert Species("H+").basename == "H"
    assert Species("N2D+").basename == "N2D"


def test_binding_energy():
    Species.reset()
    assert Species("#H").binding_energy == 600.0
    assert Species("#CH4").binding_energy == 1090.0
    assert Species("#CH4").eb == 1090.0  # test alias


def test_charge():
    Species.reset()
    assert Species("H2").charge == 0
    assert Species("H+").charge == 1
    assert Species("H-").charge == -1
    assert Species("Si++").charge == 2
    assert Species("Si++++").charge == 4


def test_atom():
    Species.reset()
    Species.set_dust_species(["GRAIN0", "GRAIN-"])
    assert Species("H").is_atom
    assert Species("C").is_atom
    assert Species("GRAIN0").is_atom
    assert not Species("E").is_atom
    assert not Species("CO").is_atom
    assert not Species("#H").is_atom


def test_element_counting():
    Species.reset()
    assert Species("H2").element_count.get("H") == 2
    assert Species("N2D+").element_count.get("N") == 2
    assert Species("N2D+").element_count.get("D") == 1
    assert Species("C10").element_count.get("C") == 10
    # special case
    assert Species("GRAIN0").element_count.get("GRAIN") == 1
    assert Species("GRAIN-").element_count.get("GRAIN") == 1


def test_eq():
    Species.reset()
    assert Species("E") == Species("e-")
    assert Species("E") in [Species("e-")]


def test_format():
    Species.reset()
    spec = Species("H")
    assert f"{spec}" == f"{spec.name}"
    assert f"{spec:<12}" == f"{spec.name:<12}"


def test_mass():
    Species.reset()
    assert Species("H").mass == 1.007
    assert Species("D").mass == 2.014
    assert Species("H2").mass == 1.007 * 2
    assert Species("He").mass == 4.002
    assert Species("CO").mass == 12.011 + 15.999


def test_massnumber():
    Species.reset()
    assert Species("D").massnumber == 2.0
    assert Species("D").A == 2.0  # test alias
    assert Species("H2").massnumber == 2.0
    assert Species("He").massnumber == 4.0
    assert Species("CO").massnumber == 28.0


def test_surface():
    Species.reset()
    Species.set_dust_species(["GRAIN0", "GRAIN-"])
    assert Species("#H2").is_surface == True
    assert Species("H2").is_surface == False
    # test changing surface symbol
    Species.surface_prefix = "G"
    assert Species("GH2").is_surface == True
    # special case
    assert Species("GRAIN0").is_surface == False
    Species.surface_prefix = "#"


# clean internal data after species tests finish
Species.reset()
