# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import getpass
from jinja2 import Template
import os
import pandas
import subprocess
from yaml import load

from pyiron.base.server.wrapper.sge import SunGridEngineCommands
from pyiron.base.server.wrapper.torque import TorqueCommands
from pyiron.base.server.wrapper.lsf import LsfCommands
from pyiron.base.server.wrapper.moab import MoabCommands
from pyiron.base.server.wrapper.slurm import SlurmCommands

__author__ = "Jan Janssen"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - " \
                "Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Feb 9, 2019"


class QueueAdapter(object):
    def __init__(self, directory='.'):
        self._config = self._read_config(file_name=os.path.join(directory, 'queue.yaml'))
        self._fill_queue_dict(queue_lst_dict=self._config['queues'])
        self._load_templates(queue_lst_dict=self._config['queues'], directory=directory)
        if self._config['queue_type'] == 'SGE':
            self._commands = SunGridEngineCommands()
        elif self._config['queue_type'] == 'TORQUE':
            self._commands = TorqueCommands()
        elif self._config['queue_type'] == 'SLURM':
            self._commands = SlurmCommands()
        elif self._config['queue_type'] == 'LSF':
            self._commands = LsfCommands()
        elif self._config['queue_type'] == 'MOAB':
            self._commands = MoabCommands()
        else:
            raise ValueError()

    @property
    def config(self):
        """

        Returns:
            dict:
        """
        return self._config

    @property
    def queue_list(self):
        """

        Returns:
            list:
        """
        return list(self._config['queues'].keys())

    @property
    def queue_view(self):
        """

        Returns:
            pandas.DataFrame:
        """
        return pandas.DataFrame(self._config['queues']).T.drop(['script', 'template'], axis=1)

    def submit_job(self, queue=None, job_name=None, working_directory=None, cores=None, memory_max=None,
                   run_time_max=None, command=None):
        """

        Args:
            queue (str/None):
            job_name (str/None):
            working_directory (str/None):
            cores (int/None):
            memory_max (int/None):
            run_time_max (int/None):
            command (str/None):

        Returns:
            int:
        """
        if isinstance(command, list):
            command = ''.join(command)
        queue_script = self._job_submission_template(queue=queue, job_name=job_name,
                                                     working_directory=working_directory, cores=cores,
                                                     memory_max=memory_max, run_time_max=run_time_max, command=command)
        queue_script_path = os.path.join(working_directory, 'run_queue.sh')
        with open(queue_script_path, 'w') as f:
            f.writelines(queue_script)
        out = self._execute_command(commands_lst=self._commands.submit_job_command + [queue_script_path],
                                    working_directory=working_directory, split_output=False)
        return int(out)

    def enable_reservation(self, process_id):
        """

        Args:
            process_id (int):

        Returns:
            str:
        """
        return self._execute_command(commands_lst=self._commands.enable_reservation_command + [str(process_id)],
                                     split_output=True)[0]

    def delete_job(self, process_id):
        """

        Args:
            process_id (int):

        Returns:
            str:
        """
        return self._execute_command(commands_lst=self._commands.delete_job_command + [str(process_id)],
                                     split_output=True)[0]

    def get_queue_status(self, user=None):
        """

        Args:
            user (str):

        Returns:
            pandas.DataFrame:
        """
        out = self._execute_command(commands_lst=self._commands.get_queue_status_command, split_output=False)
        df = self._commands.convert_queue_status(queue_status_output=out)
        if user is None:
            return df
        else:
            return df[df['user'] == user]

    def get_status_of_my_jobs(self):
        """

        Returns:
           pandas.DataFrame:
        """
        return self.get_queue_status(user=self._get_user())

    def get_status_of_job(self, process_id):
        """

        Args:
            process_id:

        Returns:
             str: ['running', 'pending', 'error']
        """
        df = self.get_queue_status()
        df_selected = df[df['jobid'] == process_id]['status']
        if len(df_selected) != 0:
            return df_selected.values[0]
        else:
            return None

    def check_queue_parameters(self, queue, cores=1, run_time_max=None, memory_max=None, active_queue=None):
        """

        Args:
            queue (str/None):
            cores (int):
            run_time_max (int/None):
            memory_max (int/None):
            active_queue (dict):

        Returns:
            list: [cores, run_time_max, memory_max]
        """
        if active_queue is None:
            active_queue = self._config['queues'][queue]
        cores = self._value_in_range(value=cores,
                                     value_min=active_queue['cores_min'],
                                     value_max=active_queue['cores_max'])
        run_time_max = self._value_in_range(value=run_time_max,
                                            value_max=active_queue['run_time_max'])
        memory_max = self._value_in_range(value=memory_max,
                                          value_max=active_queue['memory_max'])
        return cores, run_time_max, memory_max

    def _job_submission_template(self, queue=None, job_name=None, working_directory=None, cores=None, memory_max=None,
                                 run_time_max=None, command=None):
        """

        Args:
            queue (str/None):
            job_name (str/None):
            working_directory (str/None):
            cores (int/None):
            memory_max (int/None):
            run_time_max (int/None):
            command (str/None):

        Returns:
            str:
        """
        if queue is None:
            queue = self._config['queue_primary']
        for v in [job_name, working_directory, command]:
            self._value_error_if_none(value=v)
        if queue not in self.queue_list:
            raise ValueError()
        active_queue = self._config['queues'][queue]
        cores, run_time_max, memory_max = self.check_queue_parameters(queue=None,
                                                                      cores=cores,
                                                                      run_time_max=run_time_max,
                                                                      memory_max=memory_max,
                                                                      active_queue=active_queue)
        template = active_queue['template']
        return template.render(job_name=job_name,
                               working_directory=working_directory,
                               cores=cores,
                               memory_max=memory_max,
                               run_time_max=run_time_max,
                               command=command)

    @staticmethod
    def _get_user():
        """

        Returns:
            str:
        """
        return getpass.getuser()

    @staticmethod
    def _execute_command(commands_lst, working_directory=None, split_output=True):
        """

        Args:
            commands_lst (list):
            working_directory (str):
            split_output (bool):

        Returns:
            str:
        """
        if working_directory is None:
            try:
                out = subprocess.check_output(commands_lst, stderr=subprocess.STDOUT, universal_newlines=True)
            except subprocess.CalledProcessError:
                out = None
        else:
            try:
                out = subprocess.check_output(commands_lst, cwd=working_directory, stderr=subprocess.STDOUT,
                                              universal_newlines=True)
            except subprocess.CalledProcessError:
                out = None
        if out is not None and split_output:
            return out.split('\n')
        else:
            return out

    @staticmethod
    def _read_config(file_name='queue.yaml'):
        """

        Args:
            file_name (str):

        Returns:
            dict:
        """
        with open(file_name, 'r') as f:
            return load(f)

    @staticmethod
    def _fill_queue_dict(queue_lst_dict):
        """

        Args:
            queue_lst_dict (dict):
        """
        queue_keys = ['cores_min', 'cores_max', 'run_time_max', 'memory_max']
        for queue_name, queue_dict in queue_lst_dict.items():
            for key in set(queue_keys) - set(queue_dict.keys()):
                queue_dict[key] = None

    @staticmethod
    def _load_templates(queue_lst_dict, directory='.'):
        """

        Args:
            queue_lst_dict (dict):
            directory (str):
        """
        for queue_name, queue_dict in queue_lst_dict.items():
            with open(os.path.join(directory, queue_dict['script']), 'r') as f:
                queue_dict['template'] = Template(f.read())

    @staticmethod
    def _value_error_if_none(value):
        """

        Args:
            value (str/None):
        """
        if value is None:
            raise ValueError()
        if not isinstance(value, str):
            raise TypeError()

    @staticmethod
    def _value_in_range(value, value_min=None, value_max=None):
        """

        Args:
            value (int/float/None):
            value_min (int/float/None):
            value_max (int/float/None):

        Returns:
            int/float/None:
        """
        if value is not None:
            if value_min is not None and value < value_min:
                return value_min
            if value_max is not None and value > value_max:
                return value_max
            return value
        else:
            if value_min is not None:
                return value_min
            if value_max is not None:
                return value_max
            return value
