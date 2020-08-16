# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import os
import logging
from pyiron.base.project.generic import Project
from pyiron.base.settings.generic import Settings
from pyiron.base.database.filetable import get_hamilton_from_file, get_hamilton_version_from_file, \
    get_job_status_from_file

"""
The job wrapper is called from the run_job.py script, it restores the job from hdf5 and executes it.
"""

__author__ = "Joerg Neugebauer"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"

s = Settings()


class JobWrapper(object):
    """
    The job wrapper is called from the run_job.py script, it restores the job from hdf5 and executes it.

    Args:
        working_directory (str): working directory of the job
        job_id (int/ None): job ID
        hdf5_file (str): path to the HDF5 file of the job
        h5_path (str): path inside the HDF5 file to load the job
        submit_on_remote (bool): submit to queuing system on remote host
        debug (bool): enable debug mode [True/False] (optional)
    """

    def __init__(self, working_directory, job_id=None, hdf5_file=None, h5_path=None, submit_on_remote=False,
                 debug=False, connection_string=None):
        self.working_directory = working_directory
        self._remote_flag = submit_on_remote
        if connection_string is not None:
            s.open_local_sqlite_connection(connection_string=connection_string)
        pr = Project(path=os.path.join(working_directory, '..', '..'))
        if job_id is not None:
            self.job = pr.load(int(job_id))
        else:
            projectpath = s.top_path(hdf5_file)
            if projectpath is None:
                project = os.path.dirname(hdf5_file)
            else:
                project = os.path.relpath(os.path.dirname(hdf5_file), projectpath)
            job_name = h5_path[1:]
            self.job = pr.load_from_jobpath(
                job_id=None,
                db_entry={
                    "job": job_name,
                    "subjob": h5_path,
                    "projectpath": projectpath,
                    "project": project + '/',
                    "status": get_job_status_from_file(
                        hdf5_file=hdf5_file,
                        job_name=job_name
                    ),
                    "hamilton": get_hamilton_from_file(
                        hdf5_file=hdf5_file,
                        job_name=job_name
                    ),
                    "hamversion": get_hamilton_version_from_file(
                        hdf5_file=hdf5_file,
                        job_name=job_name
                    )
                },
                convert_to_object=True)

        # setup logger
        self._logger = self.setup_logger(debug=debug)

    @staticmethod
    def setup_logger(debug=False):
        """
        Setup the error logger

        Args:
            debug (bool): the level of logging, enable debug mode [True/False] (optional)

        Returns:
            logger: logger object instance
        """
        logger = logging.getLogger("pyiron_log")
        logger.setLevel(logging.INFO)
        if debug:
            logger.setLevel(logging.DEBUG)
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )
        ch.setFormatter(formatter)
        logger.addHandler(ch)
        return logger

    def run(self):
        """
        The job wrapper run command, sets the job status to 'running' and executes run_if_modal().
        """
        if self._remote_flag and self.job.server.queue is not None:
            self.job.run_if_scheduler()
        else:
            self.job.run_static()


def job_wrapper_function(working_directory, job_id=None, file_path=None, submit_on_remote=False, debug=False):
    """
    Job Wrapper function - creates a JobWrapper object and calls run() on that object

    Args:
        working_directory (str): directory where the HDF5 file of the job is located
        job_id (int/ None): job id
        file_path (str): path to the HDF5 file
        debug (bool): enable debug mode
        submit_on_remote (bool): submit to queuing system on remote host
    """
    if job_id is not None:
        job = JobWrapper(
            working_directory=working_directory,
            job_id=job_id,
            submit_on_remote=submit_on_remote,
            debug=debug
        )
    elif file_path is not None:
        hdf5_file = '.'.join(file_path.split('.')[:-1]) + '.' + file_path.split('.')[-1].split('/')[0]
        h5_path = '/'.join(file_path.split('.')[-1].split('/')[1:])
        job = JobWrapper(
            working_directory,
            job_id=None,
            hdf5_file=hdf5_file,
            h5_path='/' + h5_path,
            submit_on_remote=submit_on_remote,
            debug=debug
        )
    else:
        raise ValueError
    job.run()
