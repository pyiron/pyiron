# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
from collections import OrderedDict
from pyiron.base.master.parallel import ParallelMaster
from pyiron.atomistics.job.atomistic import AtomisticGenericJob

__author__ = "Jan Janssen"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - " \
                "Computational Materials Design (CM) Department"
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
            raise ValueError('A structure can only be set after a reference job has been assinged.')

    def get_final_structure(self):
        return self.structure


class GenericOutput(OrderedDict):
    def __init__(self):
        super(GenericOutput, self).__init__()
