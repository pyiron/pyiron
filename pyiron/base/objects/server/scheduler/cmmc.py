# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from base.objects.server.scheduler.generic import JobScheduler
from base.objects.server.queue import Queue
from base.objects.server.shelloption import ShellOption

"""
JobScheduler class - for the CMMC cluster of the Max Planck Institute fuer Eisenforschung
"""

__author__ = "Jan Janssen"
__copyright__ = "Copyright 2017, Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


class Cmmc(JobScheduler):
    """
    Example Cluster configuration for the CMMC cluster of the Max Planck Institute fuer Eisenforschung

    Args:
        working_directory (str): working directory - None by default
        scheduler_name (str): scheduler name - None by default
        cores (int): number of cores - None by default
        runtime (int): run time in seconds - 20000 by default
        wait_for_job_id (int): wait for a specific job id - False by default
        return_que_id (bool): [True/False] - True by default
    """
    def __init__(self, working_directory=None, scheduler_name=None, cores=None, runtime=20000, wait_for_job_id=False,
                 return_que_id=True):
        self.__name__ = 'cmmc'
        queue_command = 'qsub'
        queues_available = [Queue(name='impi_hydra_small', mini_cores=1, maxi_cores=40, divisor_list=[1, 2, 4, 5, 8, 10],
                                  run_time_limit=604800),
                            Queue(name='impi_hydra', mini_cores=20, maxi_cores=4240, divisor_list=20,
                                  run_time_limit=259200),
                            Queue(name='impi_hydra_cmfe.*', mini_cores=40, maxi_cores=1280, divisor_list=40,
                                  run_time_limit=259200),
                            Queue(name='impi_hy*', mini_cores=40, maxi_cores=1280, divisor_list=40,
                                  run_time_limit=259200)]
        queue_default = 'impi_hydra_small'
        queue_script_name = 'run_queue.sh'
        queue_wrapper = ['#!/bin/bash', '', 'python run_job.py']
        queue_options = {'return_que_id': ShellOption(name='return_que_id', key='-terse'),
                         'run_time': ShellOption(name='run_time', key='-l', value_prefix='h_rt='),
                         'output_file': ShellOption(name='output_file', key='-o', value='time.out'),
                         'error_file': ShellOption(name='error_file', key='-e', value='error.out'),
                         'wait_for_prev_job': ShellOption(name='wait_for_prev_job', key='-hold_jid'),
                         'working_directory': ShellOption(name='working_directory', key='-wd'),
                         'maintain_shell': ShellOption(name='maintain_shell', key='-V'),
                         'shell': ShellOption(name='shell', key='-S', value='/bin/bash'),
                         # 'job_name': ShellOption(name='job_name', key='-N'),
                         'que_name': ShellOption(name='que_name', key='-pe'),
                         'cores': ShellOption(name='cores'),
                         'script': ShellOption(name='script')}
        super(Cmmc, self).__init__(queue_command, queues_available, queue_default, queue_script_name, queue_wrapper,
                                   queue_options)
        self.quick_config(working_directory=working_directory, scheduler_name=scheduler_name, cores=cores, runtime=runtime,
                          wait_for_job_id=wait_for_job_id, return_que_id=return_que_id)
