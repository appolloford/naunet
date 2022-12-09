import pytest
from naunet.species import Species
from naunet.reactions.reaction import Reaction
from naunet.reactiontype import ReactionType
from naunet.reactions.leedsreaction import LEEDSReaction


@pytest.fixture
def reaction(request):
    return Reaction(*request.param)


@pytest.fixture
def example_freeze_reaction():
    example = LEEDSReaction(
        "7070 C11                           "
        "GC11                                              "
        "1.00E+00     0.00     132.0    541000  7"
    )
    return example


def test_init_leedsreaction_from_file(datadir):
    with open(datadir / "rate12_HO.leeds", "r") as react_file:
        for line in react_file.readlines():
            reaction = LEEDSReaction(line)
            # check all reactions can produce a reaction rate
            assert reaction.is_empty or reaction.rateexpr(grain=None)


@pytest.mark.parametrize(
    "reaction, is_equal",
    [
        (
            (
                ["C11"],
                ["#C11"],
                5.0,
                41000.0,
                1.0,
                0.0,
                132.0,
                ReactionType.GRAIN_FREEZE,
            ),
            True,
        ),
    ],
    indirect=["reaction"],
)
def test_eq_reaction(example_freeze_reaction, reaction, is_equal):

    assert example_freeze_reaction.reactants == reaction.reactants
    assert example_freeze_reaction.products == reaction.products
    assert ReactionType.GRAIN_FREEZE == LEEDSReaction.ReactionType.LEEDS_FR

    assert (example_freeze_reaction == reaction) == is_equal
    assert (
        set([example_freeze_reaction, reaction]) == set([example_freeze_reaction])
    ) == is_equal
