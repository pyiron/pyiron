# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.
"""
CLI for various pyiron utilities.
"""

import argparse
import os
import warnings

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
    parser.add_argument(
            "-d", "--dirty", action = "store_true",
            help = "do not remove pyiron log files"
    )
    subs = parser.add_subparsers()

    parser.set_defaults(cli = lambda _: parser.error("no sub command given"))

    for name, mod in cli_modules.items():
        try:
            sub_parser = subs.add_parser(name,
                help = mod.__doc__, description = mod.__doc__,
                epilog = getattr(mod, "epilog", None),
                formatter_class = getattr(mod, "formatter",
                    argparse.HelpFormatter)
            )
            sub_parser.set_defaults(cli = mod.main)
            mod.register(sub_parser)
        except AttributeError:
            warnings.warn("module '{}' does not define main or register "
                          "function, ignoring")

    args = parser.parse_args()
    args.cli(args)

    if not args.dirty:
        os.remove("pyiron.log")
