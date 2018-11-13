# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from collections import OrderedDict
from pyiron_base.core.settings.generic import Settings
from pyiron_base.objects.generic.template import PyironObject
from pyiron_base.objects.server.runmode import Runmode
from pyiron_base.objects.server.scheduler.cmmc import Cmmc
from pyiron_base.objects.server.scheduler.localhost import Localhost
import socket
import pandas

"""
Server object class which is connected to each job containing the technical details how the job is executed.
"""

__author__ = "Jan Janssen"
__copyright__ = "Copyright 2017, Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"

s = Settings()
server_types = [Cmmc]


class Server(PyironObject):  # add the option to return the job id and the hold id to the server object
    """
    Generic Server object to handle the execution environment for the job

    Args:
        host (str): hostname of the local machine
        queue (str): queue name of the currently selected queue
        cores (int): number of cores
        run_mode (str): mode of the job ['modal', 'non_modal', 'queue', 'manual']
        new_hdf (bool): create a new HDF5 file [True/False] - default=True

    Attributes:

        .. attribute:: send_to_db

            boolean option to decide which jobs should be store in the external/public database.

        .. attribute:: structure_id

            the structure ID to be linked to an external/public database.

        .. attribute:: host

            the hostname of the current system.

        .. attribute:: queue

            the que selected for a current simulation.

        .. attribute:: cores

            the number of cores selected for the current simulation.

        .. attribute:: run_time

            the run time in seconds selected for the current simulation.

        .. attribute:: run_mode

            the run mode of the job ['modal', 'non_modal', 'queue', 'manual']

        .. attribute:: new_hdf

            defines whether a subjob should be stored in the same HDF5 file or in a new one.
    """
    def __init__(self, host=None, queue=None, cores=1, run_mode='modal', new_hdf=True):
        self._scheduler = None
        self._host = None
        self._user = s.login_user
        self._run_mode = Runmode()
        self._new_hdf = None
        self._cores = None
        self._queue_id = None
        self._send_to_db = None
        self._run_time = None
        self._accept_crash = False

        self.host = host
        self.queue = queue
        self.cores = cores
        self.run_mode = run_mode
        self.new_hdf = new_hdf

        self._structure_id = None

    @property
    def send_to_db(self):
        """
        Get the boolean option to decide which jobs should be store in the external/public database

        Returns:
            bool: [True/False]
        """
        return self._send_to_db

    @send_to_db.setter
    def send_to_db(self, send):
        """
        Set the boolean option to decide which jobs should be store in the external/public database

        Args:
            send (bool): [True/False]
        """
        self._send_to_db = send

    @property
    def accept_crash(self):
        return self._accept_crash

    @accept_crash.setter
    def accept_crash(self, accept):
        self._accept_crash = accept

    @property
    def structure_id(self):
        """
        Get the structure ID to be linked to an external/public database

        Returns:
            int: structure ID
        """
        return self._structure_id

    @structure_id.setter
    def structure_id(self, structure_id):
        """
        Set the structure ID to be linked to an external/public database

        Args:
            structure_id (int): structure ID
        """
        self._structure_id = structure_id

    @property
    def host(self):
        """
        Get the hostname of the current system.
        
        Returns:
            (str): hostname
        """
        return self._host

    @host.setter
    def host(self, new_host):
        """
        Set the hostname for the current system - this should not be done manually apart from testing purposes. 
        
        Args:
            new_host (str): hostname
        """
        if new_host:
            if isinstance(new_host, str):
                self._host = new_host
            else:
                raise TypeError('The hostname should always be a string ')
        else:
            self._host = socket.gethostname()
        self._scheduler = self._init_scheduler()

    @property
    def queue(self):
        """
        The que selected for a current simulation
        
        Returns:
            (str): schedulers_name 
        """
        if self.run_mode.queue:
            return self._scheduler.active_scheduler.__name__
        else:
            return None

    @queue.setter
    def queue(self, new_scheduler):
        """
        Set a que for the current simulation, by choosing one of the available que_names
        
        Args:
            new_scheduler (str): scheduler name
        """
        if new_scheduler:
            cores = self.cores
            self._scheduler.active_scheduler = new_scheduler
            self.run_mode.queue = True
            if self._scheduler.active_scheduler.minimum_number_of_cores < cores and \
                    cores < self._scheduler.active_scheduler.maximum_number_of_cores:
                self.cores = cores

    @property
    def queue_id(self):
        """
        Get the queue ID - the ID in the queuing system is most likely not the same as the job ID.

        Returns:
            int: queue ID
        """
        return self._queue_id

    @queue_id.setter
    def queue_id(self, qid):
        """
        Set the queue ID

        Args:
            qid (int): queue ID
        """
        self._queue_id = int(qid)

    @property
    def cores(self):
        """
        The number of cores selected for the current simulation
        
        Returns:
            (int): number of cores
        """
        if self.run_mode.queue:
            return self._scheduler.cores
        else:
            return self._cores

    @cores.setter
    def cores(self, new_cores):
        """
        The number of cores selected for the current simulation
        
        Args:
            new_cores (int): number of cores
        """
        if self.run_mode.queue:
            self._scheduler.cores = new_cores
        else:
            self._cores = new_cores

    @property
    def run_time(self):
        """
        The run time in seconds selected for the current simulation
        
        Returns:
            (int): run time in seconds
        """
        if self.run_mode.queue:
            return self._scheduler.run_time
        else:
            return self._run_time

    @run_time.setter
    def run_time(self, new_run_time):
        """
        The run time in seconds selected for the current simulation
        
        Args:
            new_run_time (int): run time in seconds
        """
        if self.run_mode.queue:
            if new_run_time:
                self._scheduler.run_time = new_run_time
        else:
            self._run_time = new_run_time

    @property
    def run_mode(self):
        """
        Get the run mode of the job
        
        Returns:
            (str): ['modal', 'non_modal', 'queue', 'manual']
        """
        return self._run_mode

    @run_mode.setter
    def run_mode(self, new_mode):
        """
        Set the run mode of the job
        
        Args:
            new_mode (str): ['modal', 'non_modal', 'queue', 'manual'] 
        """
        cores = self.cores
        self._run_mode.mode = new_mode
        self.cores = cores

    @property
    def new_hdf(self):
        """
        New_hdf5 defines whether a subjob should be stored in the same HDF5 file or in a new one.
        
        Returns:
            (bool): [True / False]

        """
        return self._new_hdf

    @new_hdf.setter
    def new_hdf(self, new_hdf_bool):
        """
        New_hdf5 defines whether a subjob should be stored in the same HDF5 file or in a new one.
        
        Args:
            new_hdf_bool (bool): [True / False]
        """
        if isinstance(new_hdf_bool, bool):
            self._new_hdf = new_hdf_bool
        else:
            raise TypeError('The new_hdf5 is a boolean property, defining whether subjobs are stored in the same file.')

    def init_scheduler_run(self, working_dir, wait_for_prev_job=None, job_id=None):
        """
        Setup the job scheduler to return the scheduler options which then can be used for the python subprocess.
        
        Args:
            working_dir (str): 
            wait_for_prev_job (int): job id to wait for
            job_id (int): 

        Returns:
            (list): list of que options, [True/ False] whether the job returns a job id or not.
        """
        self._scheduler.working_directory = working_dir
        self._scheduler.write_wrapper(job_id)
        if wait_for_prev_job:
            if self._scheduler.support_wait_for_prev_job:
                self._scheduler.wait_for_job_id = wait_for_prev_job
            else:
                raise ValueError('Waiting for the previous job to be finished is not supported with this scheduler.')
        return self._scheduler.list_scheduler_options(), self._scheduler.return_scheduler_id

    def list_queues(self):
        """
        List the available Job scheduler provided by the system.

        Returns:
            (list)
        """
        return list(self._scheduler.available_schedulers_dict().keys())

    def view_queues(self):
        """
        List the available Job scheduler provided by the system.
        
        Returns:
            (pandas.DataFrame)
        """
        que_names_lst, mini_cores_lst, max_cores_lst, run_time_lst = [], [], [], []
        for que_name, que in list(self._scheduler.available_schedulers_dict().items()): 
            que_names_lst.append(que.__name__)
            mini_cores_lst.append(que.minimum_number_of_cores)
            max_cores_lst.append(que.maximum_number_of_cores)
            run_time_lst.append(que.run_time_limit) 
        return pandas.DataFrame({'minimum cores': mini_cores_lst, 'maximum cores': max_cores_lst,
                                 'run time limit': run_time_lst}, index=que_names_lst)

    def to_hdf(self, hdf, group_name=None):
        """
        Store Server object in HDF5 file
        
        Args:
            hdf: HDF5 object
            group_name (str): node name in the HDF5 file
        """
        hdf_dict = OrderedDict()
        hdf_dict["user"] = self._user
        hdf_dict["host"] = self._host
        hdf_dict["run_mode"] = self.run_mode.mode
        hdf_dict["queue"] = self.queue
        hdf_dict["qid"] = self._queue_id
        hdf_dict["cores"] = self.cores
        hdf_dict["new_h5"] = self.new_hdf
        hdf_dict["structure_id"] = self.structure_id
        hdf_dict["run_time"] = self.run_time
        hdf_dict["accept_crash"] = self.accept_crash

        if group_name:
            with hdf.open(group_name) as hdf_group:
                hdf_group["server"] = hdf_dict
        else:
            hdf["server"] = hdf_dict

    def from_hdf(self, hdf, group_name=None):
        """
        Recover Server object in HDF5 file
        
        Args:
            hdf: HDF5 object
            group_name: node name in the HDF5 file

        """
        if group_name:
            with hdf.open(group_name) as hdf_group:
                hdf_dict = hdf_group["server"]
        else:
            hdf_dict = hdf["server"]
        self._user = hdf_dict["user"]
        self._host = hdf_dict["host"]
        self._scheduler = self._init_scheduler()
        self.run_mode = hdf_dict["run_mode"]
        if self.run_mode.queue:
            self.queue = hdf_dict["queue"]
            if "qid" in hdf_dict.keys():
                self._queue_id = hdf_dict["qid"]
            else:
                self._queue_id = None
        if "structure_id" in hdf_dict.keys():
            self._structure_id = hdf_dict["structure_id"]
        self.cores = hdf_dict["cores"]
        if "run_time" in hdf_dict.keys():
            self.run_time = hdf_dict["run_time"]
        if "accept_crash" in hdf_dict.keys():
            self.accept_crash = (hdf_dict["accept_crash"] == 1)
        self.new_hdf = (hdf_dict["new_h5"] == 1)


    def db_entry(self):
        """
        connect all the info regarding the server into a single word that can be used e.g. as entry in a database
        
        Returns:
            (str): server info as single word

        """
        if self.run_mode.queue:
            server_lst = [self._host, str(self.cores), self.queue]
        else:
            server_lst = [self._host, str(self.cores)]
        return self._user + "@" + "#".join(server_lst)

    def __del__(self):
        """
        Delete the Server object from memory
        """
        del self._scheduler
        del self._user
        del self._host
        del self._run_mode

    def _init_scheduler(self):
        """
        Internal function to initialize the Job scheduler

        Returns:
            JobScheduler object
        """
        for server in server_types:
            if server.__name__.lower() in self._host:
                return server()
        return Localhost()
