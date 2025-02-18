from __future__ import annotations
import pytest
from naunet.species import Species
from naunet.grains.grain import Grain
from naunet.reactions.reaction import Reaction
from naunet.reactiontype import ReactionType
from naunet.templateloader import NetworkInfo, TemplateLoader


@pytest.fixture
def networkinfo():
    return NetworkInfo(
        [Species("H")],
        [Species("H"), Species("H2")],
        [Reaction(["H", "H"], ["H2"], 0, 100, 1.0, 0.0, 0.0, ReactionType.GAS_TWOBODY)],
        [],
        [],
        None,
        {},
    )


def test_init_templateloader():
    tl = TemplateLoader("cvode", "sparse", "cpu")
    print(tl.templates)
