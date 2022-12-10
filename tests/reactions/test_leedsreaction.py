import pytest
from naunet.species import Species
from naunet.reactions.reaction import Reaction
from naunet.reactiontype import ReactionType
from naunet.reactions.leedsreaction import LEEDSReaction


@pytest.fixture
def example_leeds_freeze_reaction():
    rstr = "7070 C11                           GC11                                              1.00E+00     0.00     132.0    541000  7"
    return LEEDSReaction(rstr)


@pytest.fixture
def example_freeze_reaction():
    return Reaction(
        ["C11"],
        ["#C11"],
        5.0,
        41000.0,
        1.0,
        0.0,
        132.0,
        ReactionType.GRAIN_FREEZE,
    )


def test_init_leedsreaction_from_file(datadir):
    with open(datadir / "rate12_HO.leeds", "r") as react_file:
        for line in react_file.readlines():
            reaction = LEEDSReaction(line)
            # check all reactions can produce a reaction rate
            assert reaction.is_empty or reaction.rateexpr(grain=None)


@pytest.mark.parametrize(
    "ref_reaction, target_reaction, is_equal",
    [("example_leeds_freeze_reaction", "example_freeze_reaction", True)],
)
def test_eq_reaction(ref_reaction, target_reaction, is_equal, request):
    assert (
        request.getfixturevalue(ref_reaction)
        == request.getfixturevalue(target_reaction)
    ) == is_equal
    assert (
        set([request.getfixturevalue(ref_reaction)])
        == set([request.getfixturevalue(target_reaction)])
    ) == is_equal
