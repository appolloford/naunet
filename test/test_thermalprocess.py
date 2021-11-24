import naunet.thermalprocess as tp
from naunet.species import Species


def test_get_allowed_heating():
    specnamelist = ["e-", "H"]
    speclist = [Species(s) for s in specnamelist]

    tp.get_allowed_heating(speclist).keys()


def test_get_allowed_cooling():
    specnamelist = ["e-", "H"]
    speclist = [Species(s) for s in specnamelist]

    assert set(tp.get_allowed_cooling(speclist).keys()) == set(["CIC_HI", "CEC_HI"])