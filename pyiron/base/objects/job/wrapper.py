# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import logging
from pyiron.base.project import Project

"""
The job wrapper is called from the run_job.py script, it restores the job from hdf5 and executes it.
"""

__author__ = "Joerg Neugebauer"
__copyright__ = "Copyright 2017, Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


class JobWrapper(object):
    """
    The job wrapper is called from the run_job.py script, it restores the job from hdf5 and executes it.

    Args:
        working_directory (str): working directory of the job
        job_id (int): job ID
        debug (bool): enable debug mode [True/False] (optional)
    """
    def __init__(self, working_directory, job_id, debug=False):
        self.working_directory = working_directory

        pr = Project(path=working_directory)
        self.job = pr.load(job_id)

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
        logger = logging.getLogger('pyiron_log')
        logger.setLevel(logging.INFO)
        if debug:
            logger.setLevel(logging.DEBUG)
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        ch.setFormatter(formatter)
        logger.addHandler(ch)
        return logger

    def run(self):
        """
        The job wrapper run command, sets the job status to 'running' and executes run_if_modal().
        """
        self.job.run_if_modal()
