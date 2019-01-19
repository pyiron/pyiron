# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from pyiron.atomistics.structure.atoms import Atoms
from pyiron.atomistics.master.parallel import ParallelMaster
from pyiron.base.master.parallel import JobGenerator

__author__ = "Jan Janssen"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH " \
                "- Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "development"
__date__ = "Sep 1, 2018"


class StructureJobGenerator(JobGenerator):
    @property
    def parameter_list(self):
        """

        Returns:
            (list)
        """
        return list(enumerate(self._job.structure_lst))

    @staticmethod
    def job_name(parameter):
        return "struct_" + str(parameter[0])

    def modify_job(self, job, parameter):
        job.structure = parameter[1]
        return job


class StructureListMaster(ParallelMaster):
    """

    Args:
        project:
        job_name:
    """

    def __init__(self, project, job_name):
        super(StructureListMaster, self).__init__(project, job_name)
        self.__name__ = 'StructureListMaster'
        self.__version__ = '0.0.1'
        self._job_generator = StructureJobGenerator(self)
        self._structure_lst = []

    @property
    def structure_lst(self):
        return self._structure_lst

    @structure_lst.setter
    def structure_lst(self, structure_lst):
        self._structure_lst = structure_lst

    def to_hdf(self, hdf=None, group_name=None):
        """
        Store the ExampleJob object in the HDF5 File

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(StructureListMaster, self).to_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input/structures") as hdf5_input:
            for ind, struct in enumerate(self.structure_lst):
                struct.to_hdf(hdf=hdf5_input, group_name='s_' + str(ind))

    def from_hdf(self, hdf=None, group_name=None):
        """
        Restore the ExampleJob object in the HDF5 File

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(StructureListMaster, self).from_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input/structures") as hdf5_input:
            self._structure_lst = [Atoms().from_hdf(hdf5_input, group_name)
                                   for group_name in hdf5_input.list_groups()]

    def collect_output(self):
        """

        Returns:

        """
        pass
