from naunet.reactions.leedsreaction import LEEDSReaction


def test_init_leedsreaction_from_file():
    leeds_reactions = []
    with open("test/test_input/rate12_full_HO.rates", "r") as react_file:
        for line in react_file.readlines():
            leeds_reactions.append(LEEDSReaction(line))
