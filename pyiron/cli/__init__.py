"""
CLI for various pyiron utilities.
"""

import argparse

from . import ls
from . import rm
from . import install
from . import reloadfile
from . import wrapper

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

cli_modules = {
        "ls":           ls,
        "rm":           rm,
        "install":      install,
        "reloadfile":   reloadfile,
        "wrapper":      wrapper
}

def main():
    parser = argparse.ArgumentParser(prog = "pyiron", description = __doc__)
    subs = parser.add_subparsers()

    parser.set_defaults(cli = lambda _: parser.error("no sub command given"))

    for name, mod in cli_modules.items():
        mod.register(subs.add_parser(name,
            help = mod.__doc__, description = mod.__doc__
        ))

    args = parser.parse_args()
    args.cli(args)

