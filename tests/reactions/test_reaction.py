import pytest
from naunet.species import Species
from naunet.reactions.reaction import Reaction
from naunet.reactiontype import ReactionType


@pytest.fixture
def example_empty_reaction():
    return Reaction()


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
        ["He"],
        ["He+", "e-"],
        5.0,
        20.0,
        0.5,
        0.0,
        0.0,
        ReactionType.GAS_COSMICRAY,
        3,
    )


@pytest.fixture
def example_cr_reaction3():
    return Reaction(
        [Species("He")],
        [Species("He+"), Species("e-")],
        5.0,
        20.0,
        0.5,
        0.0,
        0.0,
        ReactionType.GAS_COSMICRAY,
        3,
    )


@pytest.mark.parametrize(
    "ref, reacts, prods, tmin, tmax, a, b, c, rtype, idx",
    [
        (
            "example_empty_reaction",
            [],
            [],
            -1.0,
            -1.0,
            0.0,
            0.0,
            0.0,
            ReactionType.UNKNOWN,
            -1,
        ),
        (
            "example_cr_reaction1",
            [Species("He")],
            [Species("He+"), Species("e-")],
            -9999.0,
            9999.0,
            0.0,
            0.0,
            0.0,
            ReactionType.GAS_COSMICRAY,
            -1,
        ),
        (
            "example_cr_reaction2",
            [Species("He")],
            [Species("He+"), Species("e-")],
            5.0,
            20.0,
            0.5,
            0.0,
            0.0,
            ReactionType.GAS_COSMICRAY,
            3,
        ),
        (
            "example_cr_reaction3",
            [Species("He")],
            [Species("He+"), Species("e-")],
            5.0,
            20.0,
            0.5,
            0.0,
            0.0,
            101,
            3,
        ),
    ],
)
def test_init_reaction(ref, reacts, prods, tmin, tmax, a, b, c, rtype, idx, request):
    reaction = request.getfixturevalue(ref)
    assert reaction.reactants == reacts
    assert reaction.products == prods
    assert reaction.temp_min == tmin
    assert reaction.temp_max == tmax
    assert reaction.alpha == a
    assert reaction.beta == b
    assert reaction.gamma == c
    assert reaction.reaction_type == rtype
    assert reaction.idxfromfile == idx


@pytest.mark.parametrize(
    "reaction_args, is_equal",
    [
        (
            (
                ["He", "CR"],
                ["He+", "e-"],
                -9999.0,
                9999.0,
                1.0,
                0.0,
                0.0,
                ReactionType.GAS_COSMICRAY,
            ),
            True,
        ),
        (
            (
                ["He", "CR"],
                ["He+", "e-"],
                -9999.0,
                9999.0,
                1.0,
                0.0,
                0.0,
                ReactionType.UNKNOWN,
            ),
            True,
        ),
        (
            (
                ["He", "Photon"],
                ["He+", "e-"],
                -9999.0,
                9999.0,
                1.0,
                0.0,
                0.0,
                ReactionType.GAS_COSMICRAY,
            ),
            True,  # pseudo-element is not in the reactants, not compared
        ),
        (
            (
                ["He", "CR"],
                ["He+", "e-"],
                0.0,
                9999.0,
                1.0,
                0.0,
                0.0,
                ReactionType.GAS_COSMICRAY,
            ),
            False,  # different temperature
        ),
        (
            (
                ["He", "CR"],
                ["He+", "e-"],
                -9999.0,
                9999.0,
                1.0,
                0.0,
                0.0,
                ReactionType.GAS_PHOTON,
            ),
            False,  # different reaction type
        ),
    ],
)
def test_eq_reaction(example_cr_reaction1, reaction_args, is_equal):

    reaction = Reaction(*reaction_args)

    assert (example_cr_reaction1 == reaction) == is_equal
    assert (
        set([example_cr_reaction1, reaction]) == set([example_cr_reaction1])
    ) == is_equal


def test_contains(example_cr_reaction1):

    assert Species("He") in example_cr_reaction1
    assert Species("He+") in example_cr_reaction1
    assert Species("e") in example_cr_reaction1
    assert Species("H") not in example_cr_reaction1


def test_reaction_str_format(example_cr_reaction1):

    longstr = (
        "He               -> He+ + e-                        "
        ", -9999.0 < T <  9999.0, Type: GAS_COSMICRAY            "
        ", Rate: 0.0 * zeta                                                  "
        ", Source: unknown, Index: -1"
    )
    minimalstr = "He -> He+ + e-"
    shortstr = "He -> He+ + e-, -9999.0 < T <  9999.0, Type: GAS_COSMICRAY"
    naunetstr = (
        "-1   ,          He,            ,            ,"
        "         He+,          e-,            ,            ,            ,"
        " 0.000e+00, 0.000e+00, 0.000e+00,"
        " -9999.00,  9999.00, 101, unknown"
    )
    kidastr = (
        "He                                "
        "He+        e-                                            "
        "0.000e+00  0.000e+00  0.000e+00 xxxxxxxx xxxxxxxx xxxx  x  "
        "-9999   9999 xx    -1 x  x"
    )
    kromestr = "-1,He,,,He+,e-,,,,-9999.00,9999.00"
    uclchemstr = "He,CRP,NAN,He+,e-,NAN,NAN,0.0e+00,0.0,0.0,-9999.0,9999.0"

    assert str(example_cr_reaction1) == longstr
    assert f"{example_cr_reaction1}" == longstr
    assert f"{example_cr_reaction1:minimal}" == minimalstr
    assert f"{example_cr_reaction1:short}" == shortstr
    assert f"{example_cr_reaction1:naunet}" == naunetstr
    assert f"{example_cr_reaction1:kida}" == kidastr
    assert f"{example_cr_reaction1:krome}" == kromestr
    assert f"{example_cr_reaction1:uclchem}" == uclchemstr


def test_rateexpr(example_cr_reaction1):
    rate = "0.0 * zeta"
    assert example_cr_reaction1.rateexpr() == rate
