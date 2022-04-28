from naunet.reactions.leedsreaction import LEEDSReaction


def test_init_leedsreaction_from_file(datadir):
    with open(datadir / "rate12_HO.leeds", "r") as react_file:
        for line in react_file.readlines():
            reaction = LEEDSReaction(line)
            rate = reaction.rateexpr() if not reaction.is_empty else "0.0"
            assert rate  # check all reactions can produce a reaction rate
