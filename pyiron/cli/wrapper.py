"""
Run a job from hdf5.
"""

import sys
import getopt
from pyiron.base.job.wrapper import job_wrapper_function

def register(parser):
    parser.add_argument(
            "-d", "--debug", action = "store_true",
            help = "enable debug mode" # TODO: what's that mean?
    )
    parser.add_argument(
            "-j", "--job-id",
            help = "job id to run"
    )
    parser.add_argument(
            "-p", "--project",
            help = "directory where the HDF5 file of the job is located"
    )
    parser.add_argument(
            "-f", "--file-path",
            help = "path to the HDF5 file"
    )
    parser.add_argument(
            "-s", "--submit", action = "store_true",
            help = "submit to queuing system on remote host"
    )
    parser.set_defaults(cli=main)

def main(args):
    job_wrapper_function(
        working_directory=args.project,
        job_id=args.job_id,
        file_path=args.file_path,
        debug=args.debug,
        submit_on_remote=args.submit
    )
