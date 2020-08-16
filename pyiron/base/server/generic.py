# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from collections import OrderedDict
from pyiron.base.settings.generic import Settings
from pyiron.base.generic.template import PyironObject
from pyiron.base.server.runmode import Runmode
import socket

"""
Server object class which is connected to each job containing the technical details how the job is executed.
"""

__author__ = "Jan Janssen"
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


class Server(
    PyironObject
):  # add the option to return the job id and the hold id to the server object
    """
    Generic Server object to handle the execution environment for the job

    Args:
        host (str): hostname of the local machine
        queue (str): queue name of the currently selected queue
        cores (int): number of cores
        run_mode (pyiron.base.server.runmode.Runmode): mode of the job ['modal', 'non_modal', 'queue', 'manual']
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

    def __init__(
        self, host=None, queue=None, cores=1, threads=1, run_mode="modal", new_hdf=True
    ):
        self._cores = cores
        self._threads = threads
        self._run_time = None
        self._memory_limit = None
        self._host = self._init_host(host=host)

        self._active_queue = queue

        self._user = s.login_user
        self._run_mode = Runmode()
        self.run_mode = run_mode

        self._queue_id = None

        self._new_hdf = new_hdf
        self._send_to_db = False
        self._structure_id = None
        self._accept_crash = False

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
    def queue(self):
        """
        The que selected for a current simulation

        Returns:
            (str): schedulers_name
        """
        return self._active_queue

    @queue.setter
    def queue(self, new_scheduler):
        """
        Set a que for the current simulation, by choosing one of the available que_names

        Args:
            new_scheduler (str): scheduler name
        """
        if s.queue_adapter is not None:
            cores, run_time_max, memory_max = s.queue_adapter.check_queue_parameters(
                queue=new_scheduler,
                cores=self.cores,
                run_time_max=self.run_time,
                memory_max=self.memory_limit,
            )
            if cores != self.cores:
                self._cores = cores
                s.logger.debug(
                    "Updated the number of cores to: {}".format(cores)
                )
            if run_time_max != self.run_time:
                self._run_time = run_time_max
                s.logger.debug(
                    "Updated the run time limit to: {}".format(run_time_max)
                )
            if memory_max != self.memory_limit:
                self._memory_limit = memory_max
                s.logger.debug(
                    "Updated the memory limit to: {}".format(memory_max)
                )
            self._active_queue = new_scheduler
            self.run_mode = "queue"
        else:
            raise TypeError("No queue adapter defined.")

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
    def threads(self):
        return self._threads

    @threads.setter
    def threads(self, number_of_threads):
        self._threads = number_of_threads

    @property
    def cores(self):
        """
        The number of cores selected for the current simulation

        Returns:
            (int): number of cores
        """
        return self._cores

    @cores.setter
    def cores(self, new_cores):
        """
        The number of cores selected for the current simulation

        Args:
            new_cores (int): number of cores
        """
        if s.queue_adapter is not None and self._active_queue is not None:
            cores = s.queue_adapter.check_queue_parameters(
                queue=self.queue,
                cores=new_cores,
                run_time_max=self.run_time,
                memory_max=self.memory_limit,
            )[0]
            if cores != new_cores:
                self._cores = cores
                s.logger.debug(
                    "Updated the number of cores to: ", cores
                )
            else:
                self._cores = new_cores
        else:
            self._cores = new_cores

    @property
    def run_time(self):
        """
        The run time in seconds selected for the current simulation

        Returns:
            (int): run time in seconds
        """
        return self._run_time

    @run_time.setter
    def run_time(self, new_run_time):
        """
        The run time in seconds selected for the current simulation

        Args:
            new_run_time (int): run time in seconds
        """
        if s.queue_adapter is not None and self._active_queue is not None:
            run_time_max = s.queue_adapter.check_queue_parameters(
                queue=self.queue,
                cores=self.cores,
                run_time_max=new_run_time,
                memory_max=self.memory_limit,
            )[1]
            if run_time_max != new_run_time:
                self._run_time = run_time_max
                s.logger.debug(
                    "Updated the run time limit to: ", run_time_max
                )
            else:
                self._run_time = new_run_time
        else:
            self._run_time = new_run_time

    @property
    def memory_limit(self):
        return self._memory_limit

    @memory_limit.setter
    def memory_limit(self, limit):
        if s.queue_adapter is not None and self._active_queue is not None:
            memory_max = s.queue_adapter.check_queue_parameters(
                queue=self.queue,
                cores=self.cores,
                run_time_max=self.run_time,
                memory_max=limit,
            )[2]
            if memory_max != limit:
                self._memory_limit = memory_max
                s.logger.debug(
                    "Updated the memory limit to: ", memory_max
                )
            else:
                self._memory_limit = limit
        else:
            self._memory_limit = limit

    @property
    def run_mode(self):
        """
        Get the run mode of the job

        Returns:
            (str/pyiron.base.server.runmode.Runmode): ['modal', 'non_modal', 'queue', 'manual']
        """
        return self._run_mode

    @run_mode.setter
    def run_mode(self, new_mode):
        """
        Set the run mode of the job

        Args:
            new_mode (str): ['modal', 'non_modal', 'queue', 'manual']
        """
        self._run_mode.mode = new_mode
        if new_mode == "queue":
            if s.queue_adapter is None:
                raise TypeError("No queue adapter defined.")
            if self._active_queue is None:
                self.queue = s.queue_adapter.config["queue_primary"]

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
            raise TypeError(
                "The new_hdf5 is a boolean property, defining whether subjobs are stored in the same file."
            )

    @property
    def queue_list(self):
        """
        List the available Job scheduler provided by the system.

        Returns:
            (list)
        """
        return self.list_queues()

    @property
    def queue_view(self):
        """
        List the available Job scheduler provided by the system.

        Returns:
            (pandas.DataFrame)
        """
        return self.view_queues()

    @staticmethod
    def list_queues():
        """
        List the available Job scheduler provided by the system.

        Returns:
            (list)
        """
        if s.queue_adapter is not None:
            return s.queue_adapter.queue_list
        else:
            return None

    @staticmethod
    def view_queues():
        """
        List the available Job scheduler provided by the system.

        Returns:
            (pandas.DataFrame)
        """
        if s.queue_adapter is not None:
            return s.queue_adapter.queue_view
        else:
            return None

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
        hdf_dict["threads"] = self.threads
        hdf_dict["new_h5"] = self.new_hdf
        hdf_dict["structure_id"] = self.structure_id
        hdf_dict["run_time"] = self.run_time
        hdf_dict["memory_limit"] = self.memory_limit
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
        self._run_mode.mode = hdf_dict["run_mode"]
        if self.run_mode.queue:
            self._active_queue = hdf_dict["queue"]
            if "qid" in hdf_dict.keys():
                self._queue_id = hdf_dict["qid"]
            else:
                self._queue_id = None
        if "structure_id" in hdf_dict.keys():
            self._structure_id = hdf_dict["structure_id"]
        self._cores = hdf_dict["cores"]
        if "run_time" in hdf_dict.keys():
            self._run_time = hdf_dict["run_time"]
        if "memory_limit" in hdf_dict.keys():
            self._memory_limit = hdf_dict["memory_limit"]
        if "accept_crash" in hdf_dict.keys():
            self._accept_crash = hdf_dict["accept_crash"] == 1
        if "threads" in hdf_dict.keys():
            self._threads = hdf_dict["threads"]
        self._new_hdf = hdf_dict["new_h5"] == 1

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
        del self._cores
        del self._threads
        del self._run_time
        del self._memory_limit
        del self._host
        del self._active_queue
        del self._user
        del self._run_mode
        del self._queue_id
        del self._new_hdf
        del self._send_to_db
        del self._structure_id
        del self._accept_crash

    @staticmethod
    def _init_host(host):
        if host is None:
            return socket.gethostname()
        else:
            return host
