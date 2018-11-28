# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from pyiron.base.server.scheduler.generic import JobScheduler
from pyiron.base.server.queue import Queue
from pyiron.base.server.shelloption import ShellOption

"""
JobScheduler class - for any localhost
"""

__author__ = "Jan Janssen"
__copyright__ = "Copyright 2017, Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


class Localhost(JobScheduler):
    """
    Example Cluster configuration for an local host system without a real queuing system, mainly designed for
    development and testing purposes.

    Args:
        working_directory (str): working directory - None by default
        scheduler_name (str): scheduler name - None by default
        cores (int): number of cores - None by default
        runtime (int): run time in seconds - None by default
        wait_for_job_id (int): wait for a specific job id - None by default
        return_que_id (bool): [True/False] - None by default
    """
    def __init__(self, working_directory=None, scheduler_name=None, cores=None, runtime=None, wait_for_job_id=None,
                 return_que_id=None):
        queue_command = 'python'
        queues_available = [Queue(name='default', mini_cores=1, maxi_cores=1, divisor_list=1, run_time_limit=500000)]
        queue_default = 'default'
        queue_script_name = 'run_job.py'
        queue_wrapper = []
        queue_options = {'script': ShellOption(name='script')}
        super(Localhost, self).__init__(queue_command, queues_available, queue_default, queue_script_name,
                                        queue_wrapper, queue_options)
        self.quick_config(working_directory=working_directory, scheduler_name=scheduler_name, cores=cores,
                          runtime=runtime, wait_for_job_id=wait_for_job_id, return_que_id=return_que_id)
