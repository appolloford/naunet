import pytest
from naunet.species import Species
from naunet.reactions.reaction import Reaction
from naunet.reactiontype import ReactionType
from naunet.reactions.kidareaction import KIDAReaction


@pytest.fixture
def example_cr_reaction1():
    return Reaction(
        ["He", "CR"],
        ["He+", "e-"],
        -9999.0,
        9999.0,
        reaction_type=ReactionType.GAS_COSMICRAY,
    )


@pytest.fixture
def example_cr_reaction2():
    return Reaction(
        ["He", "CR"],
        ["He+", "e-"],
        -9999.0,
        9999.0,
        reaction_type=ReactionType.GAS_COSMICRAY,
        idxfromfile=999,
    )


@pytest.fixture
def example_cr_reaction3():
    return Reaction(
        ["He", "CR"],
        ["He+", "e-"],
        0.0,
        9999.0,
        reaction_type=ReactionType.GAS_COSMICRAY,
    )


@pytest.fixture
def example_kida_cr_reaction():
    rstr = "He         CR                     He+        e-                                            5.000e-01  0.000e+00  0.000e+00 2.00e+00 0.00e+00 logn  1  -9999   9999  1     3 1  1 "
    return KIDAReaction(rstr)


def test_init_kidareaction(example_kida_cr_reaction):

    reaction = example_kida_cr_reaction
    assert reaction.reactants == [Species("He")]
    assert reaction.products == [Species("He+"), Species("e-")]
    assert reaction.temp_min == -9999.0
    assert reaction.temp_max == 9999.0
    assert reaction.alpha == 0.5
    assert reaction.beta == 0.0
    assert reaction.gamma == 0.0
    assert reaction.reaction_type == ReactionType.GAS_COSMICRAY
    assert reaction.source == "kida"
    assert reaction.idxfromfile == 3


@pytest.mark.parametrize(
    "ref_reaction, target_reaction, is_equal",
    [
        ("example_kida_cr_reaction", "example_cr_reaction1", True),
        ("example_kida_cr_reaction", "example_cr_reaction2", True),
        ("example_kida_cr_reaction", "example_cr_reaction3", False),
    ],
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
    assert (
        set(
            [
                request.getfixturevalue(ref_reaction),
                request.getfixturevalue(target_reaction),
            ]
        )
        == set([request.getfixturevalue(target_reaction)])
    ) == is_equal


@pytest.mark.slow
def test_init_kidareaction_from_file(datadir):
    with open(datadir / "deuspin.kida", "r") as react_file:
        for line in react_file.readlines():
            reaction = KIDAReaction(line)
            # check all non-empty reactions produce a reaction rate
            assert reaction.is_empty or reaction.rateexpr()
