# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
from collections import OrderedDict
import inspect
from pyiron.atomistics.structure.atoms import Atoms
from pyiron.base.master.parallel import ParallelMaster, JobGenerator
from pyiron.atomistics.job.atomistic import AtomisticGenericJob
from pyiron.base.master.generic import get_function_from_string

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


class AtomisticParallelMaster(ParallelMaster, AtomisticGenericJob):
    def __init__(self, project, job_name):
        super(AtomisticParallelMaster, self).__init__(project, job_name=job_name)

    @property
    def structure(self):
        if self.ref_job:
            return self._ref_job.structure
        else:
            return None

    @structure.setter
    def structure(self, basis):
        if self.ref_job:
            self._ref_job.structure = basis
        else:
            raise ValueError(
                "A structure can only be set after a reference job has been assinged."
            )

    def get_structure(self, iteration_step=-1):
        if iteration_step == 0:
            return self.structure
        else:
            raise ValueError("iteration_step should be either 0.")


class GenericOutput(OrderedDict):
    def __init__(self):
        super(GenericOutput, self).__init__()


class MapMaster(AtomisticParallelMaster):
    def __init__(self, project, job_name):
        """

        Args:
            project:
            job_name:
        """
        super(MapMaster, self).__init__(project, job_name)
        self.__name__ = "MapMaster"
        self.__version__ = "0.0.1"
        self._job_generator = MapJobGenerator(self)
        self._map_function = None
        self.parameter_list = []

    @property
    def modify_function(self):
        return self._map_function

    @modify_function.setter
    def modify_function(self, funct):
        self._map_function = funct

    def to_hdf(self, hdf=None, group_name=None):
        """
        Store the ParameterMaster in an HDF5 file

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(MapMaster, self).to_hdf(hdf=hdf, group_name=group_name)
        if len(self.parameter_list) != 0:
            with self.project_hdf5.open("input") as hdf5_input:
                first_element = self.parameter_list[0]
                if isinstance(first_element, Atoms):
                    with hdf5_input.open("structures") as hdf5_input_str:
                        for ind, struct in enumerate(self.parameter_list):
                            struct.to_hdf(
                                hdf=hdf5_input_str, group_name="s_" + str(ind)
                            )
                elif isinstance(first_element, (int, float, str, list)):
                    hdf5_input["parameters_list"] = self.parameter_list
                else:
                    raise TypeError()
                if self._map_function is not None:
                    try:
                        hdf5_input["map_function"] = inspect.getsource(
                            self._map_function
                        )
                    except IOError:
                        hdf5_input["map_function"] = "None"
                else:
                    hdf5_input["map_function"] = "None"

    def from_hdf(self, hdf=None, group_name=None):
        """
        Restore the ParameterMaster from an HDF5 file

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(MapMaster, self).from_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf5_input:
            if "structures" in hdf5_input.list_groups():
                with hdf5_input.open("structures") as hdf5_input_str:
                    self.parameter_list = [
                        Atoms().from_hdf(hdf5_input_str, group_name)
                        for group_name in sorted(hdf5_input_str.list_groups())
                    ]
            else:
                self.parameter_list = hdf5_input["parameters_list"]
            function_str = hdf5_input["map_function"]
            if function_str == "None":
                self._map_function = None
            else:
                self._map_function = get_function_from_string(function_str)

    def collect_output(self):
        pass


class MapJobGenerator(JobGenerator):
    @property
    def parameter_list(self):
        """

        Returns:
            (list)
        """
        return self._job.parameter_list

    def modify_job(self, job, parameter):
        return self._job._map_function(job, parameter)


def pipe(project, job, step_lst, delete_existing_job=False):
    """
    Create a job pipeline

    Args:
        project (pyiron.Project): Project to calculate pipeline in
        job (AtomisticGenericJob): Template for the calculation
        step_lst (list): List of functions which create calculations

    Returns:
        FlexibleMaster:
    """
    job_lst_master = project.create_job(project.job_type.FlexibleMaster, 'lstmaster', delete_existing_job=delete_existing_job)
    if job_lst_master.status.finished:
        return job_lst_master
    else:
        for i, step_funct in enumerate(step_lst):
            job_lst_master.append(step_funct(job))
            if i > 0 and 'for_each_structure' in step_funct.__name__:
                job_lst_master.function_lst.append(_structure_many_to_many)
            elif i > 0:
                job_lst_master.function_lst.append(_structure_one_to_one)
        return job_lst_master


def _structure_one_to_one(job_prev, job_next):
    """
    Copy final structure from previous job to the input structure of the next job

    Args:
        job_prev (AtomisticGenericJob): previous Job object
        job_next (AtomisticGenericJob): next Job object
    """
    job_next.structure = job_prev.get_structure(iteration_step=-1)


def _structure_many_to_many(job_prev, job_next):
    """
    Copy list of final structures from previous job to the input structure list of the next job

    Args:
        job_prev (AtomisticParallelMaster): previous Job object
        job_next (AtomisticParallelMaster): next Job object
    """
    job_next.structure_lst = job_prev.list_structures()
