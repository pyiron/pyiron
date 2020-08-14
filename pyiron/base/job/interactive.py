# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
from pyiron.base.job.generic import GenericJob
import warnings

"""
InteractiveBase class extends the Generic Job class with all the functionality to run the job object interactivley.
"""

__author__ = "Osamu Waseda, Jan Janssen"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2018"


class InteractiveBase(GenericJob):
    """
    InteractiveBase class extends the Generic Job class with all the functionality to run the job object interactively.
    From this class all interactive Hamiltonians are derived. Therefore it should contain the properties/routines common
    to all interactive jobs. The functions in this module should be as generic as possible.

    Args:
        project (ProjectHDFio): ProjectHDFio instance which points to the HDF5 file the job is stored in
        job_name (str): name of the job, which has to be unique within the project

    Attributes:

        .. attribute:: job_name

            name of the job, which has to be unique within the project

        .. attribute:: status

            execution status of the job, can be one of the following [initialized, appended, created, submitted, running,
                                                                      aborted, collect, suspended, refresh, busy, finished]

        .. attribute:: job_id

            unique id to identify the job in the pyiron database

        .. attribute:: parent_id

            job id of the predecessor job - the job which was executed before the current one in the current job series

        .. attribute:: master_id

            job id of the master job - a meta job which groups a series of jobs, which are executed either in parallel or in
            serial.

        .. attribute:: child_ids

            list of child job ids - only meta jobs have child jobs - jobs which list the meta job as their master

        .. attribute:: project

            Project instance the jobs is located in

        .. attribute:: project_hdf5

            ProjectHDFio instance which points to the HDF5 file the job is stored in

        .. attribute:: job_info_str

            short string to describe the job by it is job_name and job ID - mainly used for logging

        .. attribute:: working_directory

            working directory of the job is executed in - outside the HDF5 file

        .. attribute:: path

            path to the job as a combination of absolute file system path and path within the HDF5 file.

        .. attribute:: version

            Version of the hamiltonian, which is also the version of the executable unless a custom executable is used.

        .. attribute:: executable

            Executable used to run the job - usually the path to an external executable.

        .. attribute:: library_activated

            For job types which offer a Python library pyiron can use the python library instead of an external executable.

        .. attribute:: server

            Server object to handle the execution environment for the job.

        .. attribute:: queue_id

            the ID returned from the queuing system - it is most likely not the same as the job ID.

        .. attribute:: logger

            logger object to monitor the external execution and internal pyiron warnings.

        .. attribute:: restart_file_list

            list of files which are used to restart the calculation from these files.

        .. attribute:: job_type

            Job type object with all the available job types: ['ExampleJob', 'SerialMaster', 'ParallelMaster', 'ScriptJob',
                                                               'ListMaster']
    """

    def __init__(self, project, job_name):
        super(InteractiveBase, self).__init__(project, job_name)
        self._interactive_library = None
        self._interactive_write_input_files = False
        self._interactive_flush_frequency = 10000
        self._interactive_write_frequency = 1
        self.interactive_cache = {}

    @property
    def interactive_flush_frequency(self):
        return self._interactive_flush_frequency

    @interactive_flush_frequency.setter
    def interactive_flush_frequency(self, frequency):
        if not isinstance(frequency, int) or frequency < 1:
            raise AssertionError("interactive_flush_frequency must be an integer>0")
        if frequency < self._interactive_write_frequency:
            raise ValueError(
                "interactive_flush_frequency must be larger or equal to interactive_write_frequency"
            )
        self._interactive_flush_frequency = frequency

    @property
    def interactive_write_frequency(self):
        return self._interactive_write_frequency

    @interactive_write_frequency.setter
    def interactive_write_frequency(self, frequency):
        if not isinstance(frequency, int) or frequency < 1:
            raise AssertionError("interactive_write_frequency must be an integer>0")
        if self._interactive_flush_frequency < frequency:
            self.interactive_flush_frequency = frequency
        self._interactive_write_frequency = frequency

    def validate_ready_to_run(self):
        """
        This should work but doesn't...
        """
        if self._interactive_flush_frequency < self._interactive_write_frequency:
            raise ValueError(
                "interactive_write_frequency must be smaller or equal to interactive_flush_frequency"
            )

    def _run_if_running(self):
        """

        Returns:

        """
        if self.server.run_mode.interactive:
            self.run_if_interactive()
        elif self.server.run_mode.interactive_non_modal:
            self.run_if_interactive_non_modal()
        else:
            super(InteractiveBase, self)._run_if_running()

    def _check_if_input_should_be_written(self):
        return (
            super(InteractiveBase, self)._check_if_input_should_be_written()
            or self._interactive_write_input_files
        )

    def interactive_is_activated(self):
        """

        Returns:

        """
        if self._interactive_library is None:
            return False
        else:
            return True

    @staticmethod
    def _extend_hdf(h5, path, key, data):
        """

        Args:
            h5:
            path:
            key:
            data:

        Returns:

        """
        if path in h5.list_groups() and key in h5[path].list_nodes():
            current_hdf = h5[path + "/" + key]
            if isinstance(data, list):
                entry = current_hdf.tolist() + data
            else:
                entry = current_hdf.tolist() + data.tolist()
            data = np.array(entry)
        h5[path + "/" + key] = data

    @staticmethod
    def _include_last_step(array, step=1, include_last=False):
        """

        Args:
            array:
            step:
            include_last:

        Returns:

        """
        if step == 1:
            return array
        if len(array) > 0:
            if len(array) > step:
                new_array = array[::step]
                index_lst = list(range(len(array)))
                if include_last and index_lst[-1] != index_lst[::step][-1]:
                    new_array.append(array[-1])
                return new_array
            else:
                if include_last:
                    return [array[-1]]
                else:
                    return []
        return []

    def interactive_flush(self, path="interactive", include_last_step=False):
        """

        Args:
            path:
            include_last_step:

        Returns:

        """
        with self.project_hdf5.open("output") as h5:
            for key in self.interactive_cache.keys():
                if len(self.interactive_cache[key]) == 0:
                    continue
                data = self._include_last_step(
                    array=self.interactive_cache[key],
                    step=self.interactive_write_frequency,
                    include_last=include_last_step,
                )
                if (
                    len(data) > 0
                    and isinstance(data[0], list)
                    and len(np.shape(data)) == 1
                ):
                    self._extend_hdf(h5=h5, path=path, key=key, data=data)
                elif np.array(data).dtype == np.dtype("O"):
                    self._extend_hdf(h5=h5, path=path, key=key, data=data)
                else:
                    self._extend_hdf(h5=h5, path=path, key=key, data=np.array(data))
                self.interactive_cache[key] = []

    def interactive_open(self):
        """

        Returns:

        """
        self.server.run_mode.interactive = True

    def interactive_close(self):
        """

        Returns:

        """
        if (
            len(list(self.interactive_cache.keys())) > 0
            and len(self.interactive_cache[list(self.interactive_cache.keys())[0]]) != 0
        ):
            self.interactive_flush(path="interactive", include_last_step=True)
        self.project_hdf5.rewrite_hdf5(job_name=self.job_name, exclude_groups=[])
        self.project.db.item_update(self._runtime(), self._job_id)
        self.status.finished = True
        self._interactive_library = None
        for key in self.interactive_cache.keys():
            self.interactive_cache[key] = []

    def interactive_store_in_cache(self, key, value):
        """

        Args:
            key:
            value:

        Returns:

        """
        self.interactive_cache[key] = value

    # def __del__(self):
    #     self.interactive_close()

    def run_if_interactive(self):
        raise NotImplementedError("run_if_interactive() is not implemented!")

    def run_if_interactive_non_modal(self):
        raise NotImplementedError("run_if_interactive_non_modal() is not implemented!")

    def to_hdf(self, hdf=None, group_name=None):
        """
        Store the InteractiveBase object in the HDF5 File

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(InteractiveBase, self).to_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf5_input:
            hdf5_input["interactive"] = {
                "interactive_flush_frequency": self._interactive_flush_frequency,
                "interactive_write_frequency": self._interactive_write_frequency,
            }

    def from_hdf(self, hdf=None, group_name=None):
        """
        Restore the InteractiveBase object in the HDF5 File

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(InteractiveBase, self).from_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf5_input:
            if "interactive" in hdf5_input.list_nodes():
                interactive_dict = hdf5_input["interactive"]
                self._interactive_flush_frequency = interactive_dict[
                    "interactive_flush_frequency"
                ]
                if "interactive_write_frequency" in interactive_dict.keys():
                    self._interactive_write_frequency = interactive_dict[
                        "interactive_write_frequency"
                    ]
                else:
                    self._interactive_write_frequency = 1
