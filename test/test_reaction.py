import pathlib

import pytest
from naunet.reactions.kidareaction import KIDAReaction
from naunet.reactions.leedsreaction import LEEDSReaction
from naunet.reactions.umistreaction import UMISTReaction


@pytest.mark.slow
def test_init_kidareaction():
    kida_reactions = []
    with open("test/test_input/deuspin.kida.uva.2017.in", "r") as react_file:
        for line in react_file.readlines():
            kida_reactions.append(KIDAReaction(line))


def test_init_leedsreaction():
    leeds_reactions = []
    with open("test/test_input/rate12_full_HO.rates", "r") as react_file:
        for line in react_file.readlines():
            leeds_reactions.append(LEEDSReaction(line))
    LEEDSReaction.finalize()


def test_init_uminstreaction():
    umist_reaction = []
    with open("test/test_input/rate12.rates", "r") as react_file:
        for line in react_file.readlines():
            umist_reaction.append(UMISTReaction(line))
