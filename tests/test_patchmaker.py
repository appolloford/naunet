from __future__ import annotations
import pytest
from naunet.species import Species
from naunet.grains.grain import Grain
from naunet.reactions.reaction import Reaction
from naunet.reactions.reactiontype import ReactionType
from naunet.templateloader import NetworkInfo
from naunet.patches import EnzoPatch


@pytest.fixture
def networkinfo():
    return NetworkInfo(
        [Species("C")],
        [Species("C"), Species("C2")],
        [Reaction(["C", "C"], ["C2"], 0, 100, 1.0, 0.0, 0.0, ReactionType.GAS_TWOBODY)],
        [],
        [],
        None,
        {},
        {},
        {},
        [],
    )


def test_init_enzopatch(networkinfo):
    patch = EnzoPatch(networkinfo, "cpu")
    print(patch.templates)
    # patch.render(["hydro_rk/Grid_CollapseMHD3DInitializeGrid.C.j2"], path="./")
