# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from .command import Command
from ...species import Species
from ...configuration import Configuration


class ExtendCommand(Command):
    """
    Extend an existing chemical network

    extend
        {input : The input filename}
        {output : The output filename}
        {--input-format=naunet : The format of the input file}
        {--output-format=naunet : The format of the output file}
        {--remove-duplicate : Remove the duplicate reactions from the network}
    """

    def __init__(self):
        super(ExtendCommand, self).__init__()

    def handle(self):

        inp = self.argument("input")
        out = self.argument("output")
        informat = self.option("input-format")
        outformat = self.option("output-format")

        from naunet.network import Network, supported_reaction_class

        for fmt in [informat, outformat]:
            if not supported_reaction_class.get(fmt):
                from importlib import util

                spec = util.spec_from_file_location(fmt, f"{fmt}.py")
                module = util.module_from_spec(spec)
                spec.loader.exec_module(module)

        net = Network(filelist=inp, fileformats=informat)

        if self.option("remove-duplicate"):
            dupes = net.find_duplicate_reaction()
            net.remove_reaction([dup[0] for dup in dupes])

        net.reindex()
        net.write(out, outformat)
