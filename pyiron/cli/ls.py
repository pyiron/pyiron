# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.
"""
Filter and prints jobs from pyiron project.
"""

import argparse
import datetime
import re
import sys
import pandas as pd
import pyiron
import pyiron.base.job.jobstatus

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

formatter = argparse.RawDescriptionHelpFormatter
epilog =  """
examples:
    Print the run time of all finished jobs:
        pyiron ls -c job totalcputime -s finished

    Print all jobs with iron:
        pyiron ls -e Fe

    Print all jobs that successfully finished yesterday and a bit:
        pyiron ls -s finished -i 1d5h

    Print all jobs that were aborted less than 5 hours ago and match
    "spx.*restart":
        pyiron ls -n "spx.*restart" -i 5h -s aborted
"""

def register(parser):

    parser.add_argument(
            "project", default = ".", nargs = "?",
            help = "path to pyiron project"
    )

    filter = parser.add_argument_group(
            title = "filter",
            description = "select which jobs to show, all filters must be true "
                        "for a job to be listed"
    )
    filter.add_argument(
            "-r", "--recursive", action = "store_true",
            help = "recurse into sub projects"
    )
    filter.add_argument(
            "-n", "--name", default = "", type = re.compile,
            help = "job name must contain this, regex allowed"
    )
    filter.add_argument(
            "-e", "--elements", nargs = "+",
            help = "chemical elements that must be present in unit cell"
    )
    filter.add_argument(
            "-s", "--status", nargs = "+",
            choices = pyiron.base.job.jobstatus.job_status_lst,
            metavar = "status",
            help = "job status must be one of the given, one of {}".format(
                    ", ".join(pyiron.base.job.jobstatus.job_status_lst))
    )
    filter.add_argument(
            "-i", "--since",
            help = "timestop must be less then the given duration before now, "
                "must be an integer and a time unit, i.e. one of d (days), "
                "h (hours), m (minutes) and s (seconds).  Combinations are "
                "possible, see examples."
    )

    output = parser.add_argument_group(
            title = "output", description = "control output style"
    )
    columns_choices = ('id', 'status', 'chemicalformula', 'job', 'subjob',
                'project', 'projectpath', 'timestart', 'timestop',
                'totalcputime', 'computer', 'hamilton', 'hamversion',
                'parentid', 'masterid')
    output.add_argument(
            "-c", "--columns", nargs = "+",
            choices = columns_choices, metavar = "column",
            default = ["id", "status", "job", "timestart", "timestop",
                        "totalcputime"],
            help = "table columns to print, pass 'all' to print whole table, "
                "one of {}".format(", ".join(columns_choices))
    )
    output.add_argument(
            "-a", "--all", action = "store_true",
            help = "show all job attributes"
    )

def main(args):

    if args.status:
        if "status" not in args.columns:
            args.columns = args.columns + ["status"]

    if args.since:
        if "timestop" not in args.columns:
            args.columns = args.columns + ["timestop"]
        try:
            matches = re.fullmatch("(\d+d)?\w*(\d+h)?\w*(\d+m)?\w*(\d+s)?",
                                args.since).groups(default = '0x')
            since = datetime.datetime.now() - datetime.timedelta(
                    days    = int(matches[0][:-1]),
                    hours   = int(matches[1][:-1]),
                    minutes = int(matches[2][:-1]),
                    seconds = int(matches[3][:-1])
            )
        except AttributeError:
            print("ERROR: {} is not a proper time delta".format(args.since),
                file = sys.stderr)
            sys.exit(1)

    table = pyiron.Project(args.project).job_table(
        full_table = True, recursive = args.recursive,
        columns = args.columns, all_columns = args.all,
        element_lst = args.elements,
        job_name_contains = args.name
    )

    if len(table) == 0:
        sys.exit(0)

    mask = [True] * len(table)
    if args.status:
        mask &= table.loc[:, "status"].isin(args.status)
    if args.since:
        mask &= table.loc[:, "timestop"] > since

    if any(mask):
        with pd.option_context("display.expand_frame_repr", False):
            print(table[mask])
