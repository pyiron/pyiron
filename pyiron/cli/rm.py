#!/usr/bin/env python3.7
"""
Remove jobs from pyiron project or whole project.
"""

import argparse
import os
import pandas as pd
import pyiron

def main():
    parser = argparse.ArgumentParser(description = __doc__)
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

    args = parser.parse_args()

    pr = pyiron.Project(args.project)
    if args.jobs_only:
        pr.remove_jobs(recursive = args.recursive)
    else:
        pr.remove(enable = True)
        if not os.listdir(args.project):
            os.rmdir(args.project)
