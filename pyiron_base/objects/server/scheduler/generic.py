# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from collections import OrderedDict
import os
from pyiron_base.objects.server.shelloption import ShellOption

"""
JobScheduler class - template for job schedulers 
"""

__author__ = "Jan Janssen"
__copyright__ = "Copyright 2017, Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


class JobScheduler(object):
    """
    The job scheduler is an abstract class which is implemented by the individual servers and handles the
    communication between the internal job scheduling and the job scheduling provided by the server.

    Args:
        scheduler_command (str): job scheduler command
        schedulers_available (list): list of available job schedulers, each as a Que object
        scheduler_default (str): name of the default scheduler
        scheduler_script_name (str): name of the scheduler script
        scheduler_wrapper (list): content of the scheduler script
        scheduler_options (dict): scheduler command line options
    """
    def __init__(self, scheduler_command, schedulers_available, scheduler_default, scheduler_script_name,
                 scheduler_wrapper, scheduler_options):
        self._active = None
        self._working_directory = None
        self._script_name = None
        self._wrapper = None
        self._wait_for_job_id = None
        self._return_scheduler_id = False
        self._options = OrderedDict()
        self.command = scheduler_command
        self.available_lst = schedulers_available
        self.wrapper = scheduler_wrapper
        options_possible = ['return_que_id', 'run_time', 'output_file', 'error_file', 'wait_for_prev_job',
                            'working_directory', 'maintain_shell', 'shell', 'que_name', 'nodes', 'cores', 'script']

        for option in options_possible:
            if option in scheduler_options.keys():
                self._options[option] = scheduler_options[option]
            else:
                self._options[option] = ShellOption(name=option, active=False)
        self.active_scheduler = scheduler_default
        self.script_name = scheduler_script_name

    @property
    def support_wait_for_prev_job(self):
        """
        Verify that the QueSystem is able to wait for an existing job to be finished first based on the que id.

        Returns:
            bool: [True/False]
        """
        return self._feature('wait_for_prev_job')

    @property
    def support_run_time_limit(self):
        """
        Verify that the QueSystem is able to stop an existing job after it reached its runtime limit.

        Returns:
            bool: [True/False]
        """
        return self._feature('run_time')

    @property
    def support_cores_limit(self):
        """
        Verify that the QueSystem is able to execute multicore jobs.

        Returns:
            bool: [True/False]
        """
        return self._feature('cores')

    @property
    def support_working_directory(self):
        """
        Verify that the QueSystem is able to force the execution to a specific directory.

        Returns:
            bool: [True/False]
        """
        return self._feature('working_directory')

    @property
    def active_scheduler(self):
        """
        Get a que for the current simulation

        Returns:
            (Queue): Queue object
        """
        return self._active
    
    @active_scheduler.setter
    def active_scheduler(self, schedulers_name):
        """
        Set a que for the current simulation, by choosing one of the available que_names

        Args:
            schedulers_name (str): queue name
        """
        if isinstance(schedulers_name, str):
            schedulers_dict = self.available_schedulers_dict()
            if schedulers_name in schedulers_dict.keys():
                self._active = schedulers_dict[schedulers_name]
                self.cores = self._active.minimum_number_of_cores
                if self.support_run_time_limit:
                    self.run_time = self._active.run_time_limit
                self._overwrite_shell_option(shell_option_name='que_name', value=self._active.__name__)
            else:
                raise ValueError('The que is not available, choose one of the following: '
                                 + str(schedulers_dict.keys()))
        else:
            raise TypeError('The schedulers_name has to be str.')

    @property
    def cores(self):
        """
        Get the number of cores selected for the current simulation

        Returns:
            int: number of cores
        """
        if self.support_cores_limit:
            if self._check_active_scheduler():
                return self._active.cores
        else:
            return 1

    @cores.setter
    def cores(self, cores):
        """
        Set the number of cores for the current simulation

        Args:
            cores (int): number of cores - need to match the requirements of the individual que
        """
        if self.support_cores_limit:
            if self._check_active_scheduler():
                self.active_scheduler.cores = cores
                self._overwrite_shell_option(shell_option_name='cores', value=cores)
        else:
            if not cores == 1:
                raise ValueError('This scheduler only supports single settings jobs')

    @property
    def run_time(self):
        """
        Get the run_time for the current simulation

        Returns:
            int: run time in seconds
        """
        if self.support_run_time_limit:
            if self._check_active_scheduler():
                return self.active_scheduler.run_time
        else:
            return None

    @run_time.setter
    def run_time(self, new_run_time):
        """
        Set the run_time for the current simulation

        Args:
            new_run_time (int): run time in seconds
        """
        if self.support_run_time_limit:
            if self._check_active_scheduler():
                self.active_scheduler.run_time = new_run_time
                self._overwrite_shell_option(shell_option_name='run_time', value=new_run_time)
        else:
            raise ValueError('This scheduler does not support run time limits')

    @property
    def working_directory(self):
        """
        Get the working directory for the current simulation.

        Returns:
            str: absolute path of the working directory
        """
        return self._working_directory

    @working_directory.setter
    def working_directory(self, work_dir):
        """
        Set the working directory for the current simulation

        Args:
            work_dir (str): absolute path of the working directory
        """
        if isinstance(work_dir, str):
            self._working_directory = work_dir
            if self.support_working_directory:
                self._overwrite_shell_option(shell_option_name='working_directory', value=work_dir)
            self._overwrite_shell_option(shell_option_name='script', value_prefix=work_dir + '/')
        else:
            raise TypeError('The working directory should be a string.')

    @property
    def script_name(self):
        """
        Get the name of the script submitted to the queing system usually a bash script named run_queue.sh .

        Returns:
            str: script name
        """
        return self._script_name

    @script_name.setter
    def script_name(self, script_name):
        """
        Set the name of the script submitted to the queing system.

        Args:
            script_name (str): script name
        """
        if isinstance(script_name, str):
            self._script_name = script_name
            self._overwrite_shell_option(shell_option_name='script', value=script_name)
        else:
            raise TypeError('The que script name should be a string.')

    @property
    def wrapper(self):
        """
        Get the shell script which is executed by the queuing system, it usually contains a shell specification and the
        run_job.py command.

        Returns:
            list: job wrapper script splitted by lines
        """
        return self._wrapper

    @wrapper.setter
    def wrapper(self, wrapper_script):
        """
        Set the shell script which is executed by the queuing system

        Args:
            wrapper_script (str, list): job wrapper script optionally splitted by lines
        """
        if isinstance(wrapper_script, str):
            wrapper_script = wrapper_script.split('\n')
        if isinstance(wrapper_script, list):
            self._wrapper = wrapper_script
        else:
            raise TypeError('The que script should either be a list of str or a single str.')

    @property
    def wait_for_job_id(self):
        """
        Get the Queue ID the current job waits for until it is released from the queuing system

        Returns:
            int: queue ID
        """
        if self.support_wait_for_prev_job:
            return self._wait_for_job_id
        else:
            return None

    @wait_for_job_id.setter
    def wait_for_job_id(self, job_id):
        """
        Set the Queue ID the current job should waits for

        Args:
            job_id (int): queue ID
        """
        if self.support_wait_for_prev_job:
            if isinstance(job_id, bool):
                self._options['wait_for_prev_job'].active = False
            elif isinstance(job_id, int):
                self._wait_for_job_id = job_id
                self._overwrite_shell_option(shell_option_name='wait_for_prev_job', value=job_id)
            else:
                raise TypeError('The job ID to wait for should be an int or use False to disable this option.')
        else:
            raise TypeError('This scheduler does not support waiting for a previous job to be finished.')

    @property
    def return_scheduler_id(self):
        """
        Some job schedulers are able to return the internal job scheduler job ID, which can be helpful for complex
        simulation protocols.

        Returns:
            bool: [True/False]
        """
        if self._feature('return_que_id'):
            return self._return_scheduler_id
        else:
            return False

    @return_scheduler_id.setter
    def return_scheduler_id(self, return_id):
        """
        Some job schedulers are able to return the internal job scheduler job ID, which can be helpful for complex
        simulation protocols.

        Args:
            return_id (bool): [True/False]
        """
        if self._feature('return_que_id'):
            if isinstance(return_id, bool):
                self._return_scheduler_id = return_id
                self._options['return_que_id'].active = return_id
            else:
                raise TypeError('The return schedulder id is a boolean flag, it can only be true or false.')
        else:
            if return_id:
                raise TypeError('The return schedulder id is not supported for this schedulder.')

    def _feature(self, option):
        """
        internal function to check if the possible options are available for the specific job scheduler and activated.

        Args:
            option (str): one of the following ['return_que_id', 'run_time', 'output_file', 'error_file',
                         'wait_for_prev_job', 'working_directory', 'maintain_shell', 'shell', 'que_name', 'nodes',
                         'cores', 'script']

        Returns:
            bool: [True/False]
        """
        if option in self._options.keys():
            return self._options[option].active
        else:
            return False

    def _check_active_scheduler(self):
        """
        Verify that one queue is set to active for this job.

        Returns:
            bool: [True/False]
        """
        if self.active_scheduler:
            return True
        else:
            raise ValueError('First set a que, so the number of cores and runtime can be validated.')

    def _overwrite_shell_option(self, shell_option_name, value=None, value_prefix=None):
        """
        Update the list of options for the shell command

        Args:
            shell_option_name (str): name of the option
            value (str): the value of the option
            value_prefix (str): the value prefix
        """
        if shell_option_name in self._options.keys():
            prev_shell_option = self._options[shell_option_name]
            if value_prefix and value:
                self._options[shell_option_name] = \
                    ShellOption(name=prev_shell_option.__name__, key=prev_shell_option.key, value=str(value),
                                value_prefix=value_prefix, active=True)
            elif not value_prefix and value:
                self._options[shell_option_name] = \
                    ShellOption(name=prev_shell_option.__name__, key=prev_shell_option.key, value=str(value),
                                value_prefix=prev_shell_option.value_prefix, active=True)
            elif value_prefix and not value:
                self._options[shell_option_name] = \
                    ShellOption(name=prev_shell_option.__name__, key=prev_shell_option.key,
                                value=prev_shell_option.value, value_prefix=value_prefix, active=True)
            else:
                raise ValueError('Overwrite either the value or the value_prefix, or both.')

    def list_scheduler_options(self):
        """
        Generate the list of queue options including the queue command, so it can be used as input for an sub process.

        Returns:
            list: list of queue options
        """
        if self._check_active_scheduler():
            que_option_lst = [self.command]
            for name, que_option in self._options.items():
                for que_command in que_option.output_as_list():
                    que_option_lst.append(que_command)
            return que_option_lst

    def write_wrapper(self, job_id=None):
        """
        Write the wrapper script based on the information stored in self.wrapper into the script self.script_name.
        """
        if self.working_directory:
            if job_id:
                self.script_name = 'pi_'+str(job_id)+'.sh'
            file_name = os.path.join(self.working_directory, self.script_name)
            if self.wrapper:
                with open(file_name, "w") as f:
                    for line in self.wrapper:
                        f.write(line + '\n')
        else:
            raise ValueError('Working directory not set.')

    def available_schedulers_dict(self):
        """
        Display all available queues and their settings

        Returns:
            dict: dictionary of all available queues
        """
        return dict([[que.__name__, que] for que in self.available_lst])

    def quick_config(self, working_directory=None, scheduler_name=None, cores=None, runtime=None, wait_for_job_id=None,
                     return_que_id=False):
        """
        Quick configuration of the Job Scheduler, sets all necessary properties, especially useful for tests.

        Args:
            working_directory (str): working directory - None by default
            scheduler_name (str): scheduler name - None by default
            cores (int): number of cores - None by default
            runtime (int): run time in seconds - None by default
            wait_for_job_id (int): wait for a specific job id - None by default
            return_que_id (bool): [True/False] - False by default
        """
        if working_directory:
            self.working_directory = working_directory
        if scheduler_name:
            self.active_scheduler = scheduler_name
        if cores:
            self.cores = cores
        if runtime:
            self.run_time = runtime
        if wait_for_job_id:
            self.wait_for_job_id = wait_for_job_id
        else:
            if self.support_wait_for_prev_job:
                self.wait_for_job_id = False
        if return_que_id:
            self.return_scheduler_id = True
        else:
            self.return_scheduler_id = False

    def __repr__(self):
        """
        Human readable representation of the available ques and the options used by the queing system.
        """
        output_str = 'Available Ques:\n' + str(self.available_schedulers_dict()) + '\n\n' + \
                     'Que Options:\n' + str(self.list_scheduler_options())
        if self.wrapper:
            que_wrapper_str = ''
            for line in self.wrapper:
                que_wrapper_str += line + '\n'
            output_str += '\n\n Que Script:\n' + que_wrapper_str
        return output_str
