# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.
"""
Remove jobs from pyiron project or whole project.
"""

import os
import pyiron

__author__ = "Marvin Poul"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Marvin Poul"
__email__ = "poul@mpie.de"
__status__ = "production"
__date__ = "23 Jun, 2020"

def register(parser):
    parser.add_argument(
            "project", default = ".", nargs = "?",
            help = "path to pyiron project"
    )
    parser.add_argument(
            "-j", "--jobs-only", action = "store_true",
            help = "only remove jobs inside project"
    )
    parser.add_argument(
            "-r", "--recursive", action = "store_true",
            help = "recurse into subprojects"
    )

def main(args):

    pr = pyiron.Project(args.project)
    if args.jobs_only:
        pr.remove_jobs_silently(recursive = args.recursive)
    else:
        pr.remove(enable = True)
        if not os.listdir(args.project):
            os.rmdir(args.project)
