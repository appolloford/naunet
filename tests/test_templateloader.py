from __future__ import annotations
import pytest
from naunet.species import Species
from naunet.dusts.dust import Dust
from naunet.reactions.reaction import Reaction
from naunet.reactions.reactiontype import ReactionType
from naunet.templateloader import NetworkInfo, TemplateLoader


@pytest.fixture
def networkinfo():
    return NetworkInfo(
        [Species("H"), Species("H2")],
        [Reaction(["H", "H"], ["H2"], 0, 100, 1.0, 0.0, 0.0, ReactionType.GAS_TWOBODY)],
        [],
        [],
        None,
        {},
        {},
        {},
        [],
    )


def test_init_templateloader(networkinfo):
    tl = TemplateLoader(networkinfo, "cvode", "sparse", "cpu")
    print(tl.templates)
