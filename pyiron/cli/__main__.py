"""
CLI for various pyiron utilities.
"""

import argparse

from . import ls
from . import rm
from . import install

__author__ = "Marvin Poul"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Marvin Poul"
__email__ = "poul@mpie.de"
__status__ = "development"
__date__ = "26 Jun, 2020"

parser = argparse.ArgumentParser(prog = "pyiron", description = __doc__)
subs = parser.add_subparsers()

parser.set_defaults(cli = lambda _: parser.error("no sub command given"))

ls.register(subs.add_parser("ls",
        help = ls.__doc__, description = ls.__doc__
))

rm.register(subs.add_parser("rm",
        help = rm.__doc__, description = rm.__doc__
))

install.register(subs.add_parser("install",
        help = install.__doc__, description = install.__doc__
))

args = parser.parse_args()
args.cli(args)
