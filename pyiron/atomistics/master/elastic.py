# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function

import numpy as np
import scipy.constants
import warnings
from pyiron.atomistics.structure.atoms import Atoms, ase_to_pyiron
from pyiron.atomistics.master.parallel import AtomisticParallelMaster
from pyiron_base import JobGenerator

__author__ = "Joerg Neugebauer, Jan Janssen"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


eV_div_A3_to_GPa = (
    1e21 / scipy.constants.physical_constants["joule-electron volt relationship"][0]
)


class ElasticJobGenerator(JobGenerator):
    @property
    def parameter_list(self):
        """

        Returns:
            (list)
        """
        parameter_lst = []
        if len(self._job.input['strain_matrices'])==0:
            eps_lst = np.random.random((int(self._job.input['num_points']), 3, 3))-0.5
            eps_lst *= self._job.input['max_strain']
            eps_lst += np.einsum('nij->nji', eps_lst)
            self._job.input['strain_matrices'] = eps_lst.tolist()
        for ii, epsilon in enumerate(self._job.input['strain_matrices']):
            basis = self._job.ref_job.structure.copy()
            basis.apply_strain(np.array(epsilon))
            parameter_lst.append([ii, basis])
        return parameter_lst

    def job_name(self, parameter):
        return "{}_{}".format(self._job.job_name, parameter[0]).replace(".", "_")

    def modify_job(self, job, parameter):
        job.structure = parameter[1]
        return job

# ToDo: not all abstract methods implemented
class ElasticTensor(AtomisticParallelMaster):
    def __init__(self, project, job_name="elastic"):
        """

        Args:
            project:
            job_name:
        """
        super().__init__(project, job_name)
        self.__name__ = "ElasticTensor"
        self.__version__ = "0.1.0"

        # print ("h5_path: ", self.project_hdf5._h5_path)

        # define default input
        self.input["num_points"] = (11, "number of sample points")
        self.input["max_strain"] = (
            0.01,
            "relative volume variation around volume defined by ref_ham",
        )
        self.input['strain_matrices'] = []
        self._job_generator = ElasticJobGenerator(self)

    def list_structures(self):
        if self.ref_job.structure is not None:
            return [parameter[1] for parameter in self._job_generator.parameter_list]
        else:
            return []

    def collect_output(self):
        if self.ref_job.server.run_mode.interactive:
            ham = self.project_hdf5.inspect(self.child_ids[0])
            self._output["energy"] = ham["output/generic/energy_tot"]
            self._output["pressures"] = ham["output/generic/pressures"]
        else:
            erg_lst, pressure_lst, id_lst = [], [], []
            for job_id in self.child_ids:
                ham = self.project_hdf5.inspect(job_id)
                print("job_id: ", job_id, ham.status)
                if "energy_tot" in ham["output/generic"].list_nodes():
                    erg_lst.append(ham["output/generic/energy_tot"][-1])
                elif "energy_pot" in ham["output/generic"].list_nodes():
                    erg_lst.append(ham["output/generic/energy_pot"][-1])
                else:
                    raise ValueError('Neither energy_pot or energy_tot was found.')
                if "pressures" in ham['output/generic'].list_nodes():
                    pressure_lst.append(ham["output/generic/pressures"][-1])
                id_lst.append(job_id)
            self._output["pressures"] = np.array(pressure_lst)
            self._output["energy"] = np.array(erg_lst)
            self._output["id"] = np.array(id_lst)
        with self.project_hdf5.open("output") as hdf5_out:
            for key, val in self._output.items():
                hdf5_out[key] = val


