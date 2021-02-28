import pathlib

import pytest
from naunet.reactions.kidareaction import KIDAReaction
from naunet.reactions.leedsreaction import LEEDSReaction


@pytest.mark.slow
def test_init_kidareaction():
    kida_reactions = []
    with open("test/test_input/deuspin.kida.uva.2017.in", "r") as react_file:
        for line in react_file.readlines():
            kida_reactions.append(KIDAReaction(line))


def test_init_leedsreactio():
    leeds_reactions = []
    with open("test/test_input/rate12_full_HO.rates", "r") as react_file:
        for line in react_file.readlines():
            leeds_reactions.append(LEEDSReaction(line))