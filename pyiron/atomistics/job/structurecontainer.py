# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
from pyiron.atomistics.structure.atoms import Atoms
from pyiron.atomistics.job.atomistic import AtomisticGenericJob

__author__ = "Yury Lysogorskiy"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - " \
                "Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "development"
__date__ = "Sep 1, 2017"


class StructureContainer(AtomisticGenericJob):
    def __init__(self, project, job_name):
        super(StructureContainer, self).__init__(project, job_name)
        self.__version__ = "0.1"
        self.__name__ = "StructureContainer"
        self.server.run_mode.interactive = True

    @property
    def structure(self):
        return self._structure

    @structure.setter
    def structure(self, structure):
        if isinstance(structure, AtomisticGenericJob):
            if structure.structure:
                self._structure = structure.get_structure(iteration_step=-1)
                self.parent_id = structure.job_id
            else:
                raise ValueError('The job does not contain any structure to import.')
        elif isinstance(structure, Atoms):
            self._structure = structure
        else:
            raise TypeError('You can only append a structure object or an GenericJob object.')
        if self.status.initialized:
            self.run()

    def append(self, structure_to_append):
        self.structure = structure_to_append

    def run_if_interactive(self):
        self.status.finished = True

    def to_hdf(self, hdf=None, group_name=None):
        super(StructureContainer, self).to_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf5_input:
            self.structure.to_hdf(hdf5_input)

    def from_hdf(self, hdf=None, group_name=None):
        super(StructureContainer, self).from_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf5_input:
            self.structure = Atoms().from_hdf(hdf5_input)
