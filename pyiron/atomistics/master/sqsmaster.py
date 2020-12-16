# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
from pyiron.atomistics.master.parallel import AtomisticParallelMaster
from pyiron_base import JobGenerator

__author__ = "Jan Janssen"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "0.0.1"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "development"
__date__ = "Oct 29, 2020"


class SQSMaster(AtomisticParallelMaster):
    def __init__(self, project, job_name):
        super(SQSMaster, self).__init__(project, job_name)
        self.__name__ = "SQSMaster"
        self.__version__ = "0.0.1"
        self.input["fraction_lst"] = []
        self.input["species_one"] = ""
        self.input["species_two"] = ""
        self._job_generator = SQSJobGenerator(self)

    def collect_output(self):
        pass

    @property
    def list_of_structures(self):
        return [self.project_hdf5.load(job_id).list_of_structures[0] for job_id in self.child_ids]

    def list_structures(self):
        return self.list_of_structures


class SQSJobGenerator(JobGenerator):
    @property
    def parameter_list(self):
        """

        Returns:
            (list)
        """
        return [
            [
                "sqs_" + str(np.round(f, 4)).replace(".", "_"),
                {self._master.input["species_one"]: f, self._master.input["species_two"]: 1-f}
            ]
            for f in self._master.input["fraction_lst"]
        ]

    @staticmethod
    def job_name(parameter):
        return parameter[0]

    def modify_job(self, job, parameter):
        job.input['mole_fractions'] = parameter[1]
        return job
