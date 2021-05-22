from naunet.species import Species


def test_alias():
    assert Species("H2").alias == "H2I"
    assert Species("CO").alias == "COI"
    assert Species("e-").alias == "eM"
    assert Species("H-").alias == "HM"
    assert Species("H+").alias == "HII"
    assert Species("oH2D+").alias == "oH2DII"
    assert Species("N2D+").alias == "N2DII"
    assert Species("#CO").alias == "GCOI"


def test_alias_setter():
    spec = Species("H2")
    spec.alias = "HH"
    assert spec.alias == "HH"


def test_basename():
    assert Species("H2").basename == "H2"
    assert Species("CO").basename == "CO"
    assert Species("H+").basename == "H"
    assert Species("N2D+").basename == "N2D"


def test_binding_energy():
    assert Species("#H").binding_energy == 600.0
    assert Species("#CH4").binding_energy == 1090.0


def test_charge():
    assert Species("H2").charge == 0
    assert Species("H+").charge == 1
    assert Species("H-").charge == -1
    assert Species("Si++").charge == 2
    assert Species("Si++++").charge == 4


def test_element_counting():
    assert Species("H2").element_count.get("H") == 2
    assert Species("N2D+").element_count.get("N") == 2
    assert Species("N2D+").element_count.get("D") == 1
    assert Species("C10").element_count.get("C") == 10
    # special case
    assert Species("GRAIN0").element_count.get("GRAIN") == 1
    assert Species("GRAIN-").element_count.get("GRAIN") == 1


def test_mass():
    assert Species("H").mass == 1.007
    assert Species("H2").mass == 1.007 * 2
    assert Species("He").mass == 4.002
    assert Species("CO").mass == 12.011 + 15.999


def test_massnumber():
    assert Species("H2").massnumber == 2.0
    assert Species("He").massnumber == 4.0
    assert Species("CO").massnumber == 28.0


def test_surface():
    assert Species("#H2").is_surface == True
    assert Species("H2").is_surface == False
    # test changing surface symbol
    Species.surface_prefix = "G"
    assert Species("GH2").is_surface == True
    # special case
    assert Species("GRAIN0").is_surface == False