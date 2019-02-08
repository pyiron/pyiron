import getpass
from jinja2 import Template
import os
import pandas
import subprocess
from yaml import load

from pyiron.base.server.wrapper.sge import SunGridEngineCommands


class QueueAdapter(object):
    def __init__(self, directory='.'):
        self._config = self._read_config(file_name=os.path.join(directory, 'queue.yaml'))
        self._fill_queue_dict(queue_lst_dict=self._config['queues'])
        self._load_templates(queue_lst_dict=self._config['queues'], directory=directory)
        if self._config['queue_type'] == 'SGE':
            self._commands = SunGridEngineCommands()
        else:
            raise ValueError()

    @property
    def queue_list(self):
        return list(self._config['queues'].keys())

    @property
    def queue_view(self):
        return pandas.DataFrame(self._config['queues']).T.drop(['script', 'template'], axis=1)

    def submit_job(self, queue=None, job_name=None, working_directory=None, cores=None, memory_max=None,
                   run_time_max=None, command=None):
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
        return self._execute_command(commands_lst=self._commands.enable_reservation_command + [str(process_id)],
                                     split_output=True)[0]

    def delete_job(self, process_id):
        return \
        self._execute_command(commands_lst=self._commands.delete_job_command + [str(process_id)], split_output=True)[0]

    def get_queue_status(self, user=None):
        out = self._execute_command(commands_lst=self._commands.get_queue_status_command, split_output=False)
        df = self._commands.convert_queue_status(queue_status_output=out)
        if user is None:
            return df
        else:
            return df[df['user'] == user]

    def get_status_of_my_jobs(self):
        return self.get_queue_status(user=self._get_user())

    def get_status_of_job(self, process_id):
        df = self.get_queue_status()
        return df[df['jobid'] == process_id]['status'].values[0]

    def _job_submission_template(self, queue=None, job_name=None, working_directory=None, cores=None, memory_max=None,
                                 run_time_max=None, command=None):
        for v in [queue, job_name, working_directory, command]:
            self._value_error_if_none(value=v)
        if queue not in self.queue_list:
            raise ValueError()
        active_queue = self._config['queues'][queue]
        self._value_in_range(value=cores, value_min=active_queue['cores_min'], value_max=active_queue['cores_max'])
        self._value_in_range(value=run_time_max, value_max=active_queue['run_time_max'])
        self._value_in_range(value=memory_max, value_max=active_queue['memory_max'])
        template = active_queue['template']
        return template.render(job_name=job_name,
                               working_directory=working_directory,
                               cores=cores,
                               memory_max=memory_max,
                               run_time_max=run_time_max,
                               command=command)

    @staticmethod
    def _get_user():
        return getpass.getuser()

    @staticmethod
    def _execute_command(commands_lst, working_directory=None, split_output=True):
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
        with open(file_name, 'r') as f:
            return load(f)

    @staticmethod
    def _fill_queue_dict(queue_lst_dict):
        queue_keys = ['cores_min', 'cores_max', 'run_time_max', 'memory_max']
        for queue_name, queue_dict in queue_lst_dict.items():
            for key in set(queue_keys) - set(queue_dict.keys()):
                queue_dict[key] = None

    @staticmethod
    def _load_templates(queue_lst_dict, directory='.'):
        for queue_name, queue_dict in queue_lst_dict.items():
            with open(os.path.join(directory, queue_dict['script']), 'r') as f:
                queue_dict['template'] = Template(f.read())

    @staticmethod
    def _value_error_if_none(value):
        if value is None:
            raise ValueError()
        if not isinstance(value, str):
            raise TypeError()

    @staticmethod
    def _value_in_range(value, value_min=None, value_max=None):
        if value is not None:
            if value_min is not None and value < value_min:
                raise ValueError()
            if value_max is not None and value > value_max:
                raise ValueError()
