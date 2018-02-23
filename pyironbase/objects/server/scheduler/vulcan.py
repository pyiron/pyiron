# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from pyironbase.objects.server.scheduler.generic import JobScheduler
from pyironbase.objects.server.queue import Queue
from pyironbase.objects.server.shelloption import ShellOption

"""
JobScheduler class - for the Vulcan cluster of ICAMS  
"""

__author__ = "Yury Lysogorskiy"
__copyright__ = "Copyright 2017, Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Yury Lysogorskiy"
__email__ = "yury.lysogorskiy@icams.rub.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


class Vulcan(JobScheduler):
    """
    Example Cluster configuration for the Vulcan cluster of ICAMS, Ruhr-University, Bochum

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
        self.__name__ = 'vulcan'
        queue_command = 'qsub'
        queues_available = [Queue(name='smp', mini_cores=1, maxi_cores=8, divisor_list=1,
                                  run_time_limit=2 * 7 * 24 * 3600),
                            Queue(name='smp8', mini_cores=8, maxi_cores=8, divisor_list=8,
                                  run_time_limit=2 * 7 * 24 * 3600),
                            Queue(name='smp12', mini_cores=12, maxi_cores=12, divisor_list=12,
                                  run_time_limit=2 * 7 * 24 * 3600),
                            Queue(name='mpi8', mini_cores=8, maxi_cores=64, divisor_list=8,
                                  run_time_limit=2 * 7 * 24 * 3600),
                            Queue(name='mpi12', mini_cores=12, maxi_cores=256, divisor_list=12,
                                  run_time_limit=2 * 7 * 24 * 3600),
                            Queue(name='mpi16', mini_cores=16, maxi_cores=512, divisor_list=16,
                                  run_time_limit=2 * 7 * 24 * 3600),
                            Queue(name='mpi20', mini_cores=20, maxi_cores=240, divisor_list=20,
                                  run_time_limit=1 * 7 * 24 * 3600)
                            ]
        queue_default = 'smp'
        queue_script_name = 'run_queue.sh'
        queue_wrapper = ['#!/bin/bash', '', 'source ~/.bashrc', 'python run_job.py']
        queue_options = {'return_que_id': ShellOption(name='return_que_id', key='-terse'),
                         'run_time': ShellOption(name='run_time', key='-l', value_prefix='h_rt='),
                         'output_file': ShellOption(name='output_file', key='-o', value='time.out'),
                         'error_file': ShellOption(name='error_file', key='-e', value='error.out'),
                         'wait_for_prev_job': ShellOption(name='wait_for_prev_job', key='-hold_jid'),
                         'working_directory': ShellOption(name='working_directory', key='-wd'),
                         'shell': ShellOption(name='shell', key='-S', value='/bin/bash'),
                         'que_name': ShellOption(name='que_name', key='-pe'),
                         'cores': ShellOption(name='cores'),
                         'script': ShellOption(name='script')}

        super(Vulcan, self).__init__(queue_command, queues_available, queue_default, queue_script_name, queue_wrapper,
                                     queue_options)
        self.quick_config(working_directory=working_directory, scheduler_name=scheduler_name, cores=cores,
                          runtime=runtime, wait_for_job_id=wait_for_job_id, return_que_id=return_que_id)
