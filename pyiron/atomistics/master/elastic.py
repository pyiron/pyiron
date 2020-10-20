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
        for strain in np.linspace(
            1 - self._job.input["vol_range"],
            1 + self._job.input["vol_range"],
            self._job.input["num_points"],
        ):
            basis = self._job.ref_job.structure.copy()
            basis.set_cell(basis.cell * strain ** (1.0 / 3.0), scale_atoms=True)
            parameter_lst.append([np.round(strain, 7), basis])
        return parameter_lst

    @staticmethod
    def job_name(parameter):
        return "strain_" + str(parameter[0]).replace(".", "_")

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
        super(Murnaghan, self).__init__(project, job_name)
        self.__name__ = "ElasticTensor"
        self.__version__ = "0.1.0"

        # print ("h5_path: ", self.project_hdf5._h5_path)

        # define default input
        self.input["num_points"] = (11, "number of sample points")
        self.input["max_strain"] = (
            0.01,
            "relative volume variation around volume defined by ref_ham",
        )

        self._job_generator = ElasticJobGenerator(self)

    def list_structures(self):
        if self.ref_job.structure is not None:
            return [parameter[1] for parameter in self._job_generator.parameter_list]
        else:
            return []

    def collect_output(self):
        if self.ref_job.server.run_mode.interactive:
            ham = self.project_hdf5.inspect(self.child_ids[0])
            erg_lst = ham["output/generic/energy_tot"]
            vol_lst = ham["output/generic/volume"]
            arg_lst = np.argsort(vol_lst)

            self._output["volume"] = vol_lst[arg_lst]
            self._output["energy"] = erg_lst[arg_lst]
        else:
            erg_lst, vol_lst, err_lst, id_lst = [], [], [], []
            for job_id in self.child_ids:
                ham = self.project_hdf5.inspect(job_id)
                print("job_id: ", job_id, ham.status)
                if "energy_tot" in ham["output/generic"].list_nodes():
                    energy = ham["output/generic/energy_tot"][-1]
                elif "energy_pot" in ham["output/generic"].list_nodes():
                    energy = ham["output/generic/energy_pot"][-1]
                else:
                    raise ValueError('Neither energy_pot or energy_tot was found.')
                volume = ham["output/generic/volume"][-1]
                erg_lst.append(np.mean(energy))
                err_lst.append(np.var(energy))
                vol_lst.append(volume)
                id_lst.append(job_id)
            vol_lst = np.array(vol_lst)
            erg_lst = np.array(erg_lst)
            err_lst = np.array(err_lst)
            id_lst = np.array(id_lst)
            arg_lst = np.argsort(vol_lst)

            self._output["volume"] = vol_lst[arg_lst]
            self._output["energy"] = erg_lst[arg_lst]
            self._output["error"] = err_lst[arg_lst]
            self._output["id"] = id_lst[arg_lst]

        with self.project_hdf5.open("output") as hdf5_out:
            for key, val in self._output.items():
                hdf5_out[key] = val
        if self.input["fit_type"] == "polynomial":
            self.fit_polynomial(fit_order=self.input["fit_order"])
        else:
            self._fit_eos_general(fittype=self.input["fit_type"])


