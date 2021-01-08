# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
from pyiron_atomistic.atomistics.master.murnaghan import MurnaghanJobGenerator
from pyiron_atomistic.atomistics.master.parallel import AtomisticParallelMaster

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


def calc_v0_from_fit_funct(fit_funct, x, save_range=0.0, return_ind=False):
    fit_funct_der = fit_funct.deriv().r
    fit_funct_der_r = fit_funct_der[fit_funct_der.imag == 0].real
    fit_funct_der_val = fit_funct.deriv(2)(fit_funct_der_r)
    select = (fit_funct_der_val > 0) & (fit_funct_der_r > np.min(x) * (1-save_range)) & (fit_funct_der_r < np.max(x) *(1+save_range))
    v0_lst = fit_funct_der_r[select]
    if len(v0_lst) == 1:
        if return_ind:
            return v0_lst[0], (np.abs(x - v0_lst[0])).argmin()
        else:
            return v0_lst[0]
    else:
        select = fit_funct_der_val > 0
        v0_lst = fit_funct_der_r[select]
        if len(v0_lst) == 1:
            if return_ind:
                return v0_lst[0], (np.abs(x - v0_lst[0])).argmin()
            else:
                return v0_lst[0]
        else:
            if return_ind:
                return None, None
            else:
                return None


class QuasiHarmonicJob(AtomisticParallelMaster):
    def __init__(self, project, job_name="murnaghan"):
        """

        Args:
            project:
            job_name:
        """
        super(QuasiHarmonicJob, self).__init__(project, job_name)
        self.__name__ = "QuasiHarmonicJob"
        self.__version__ = "0.0.1"

        # define default input
        self.input["num_points"] = (11, "number of sample points")
        self.input["vol_range"] = (
            0.1,
            "relative volume variation around volume defined by ref_ham",
        )
        self.input["temperature_start"] = 0
        self.input["temperature_end"] = 500
        self.input["temperature_steps"] = 10
        self.input["polynomial_degree"] = 3
        self._job_generator = MurnaghanJobGenerator(self)

    def collect_output(self):
        free_energy_lst, entropy_lst, cv_lst, volume_lst = [], [], [], []
        for job_id in self.child_ids:
            job = self.project_hdf5.load(job_id)
            thermal_properties = job.get_thermal_properties(temperatures=np.linspace(
                self.input["temperature_start"],
                self.input["temperature_end"],
                int(self.input["temperature_steps"])
            ))
            free_energy_lst.append(thermal_properties.free_energies)
            entropy_lst.append(thermal_properties.entropy)
            cv_lst.append(thermal_properties.cv)
            volume_lst.append(job.structure.get_volume())

        arg_lst = np.argsort(volume_lst)

        self._output["free_energy"] = np.array(free_energy_lst)[arg_lst]
        self._output["entropy"] = np.array(entropy_lst)[arg_lst]
        self._output["cv"] = np.array(cv_lst)[arg_lst]

        temperature_mesh, volume_mesh = np.meshgrid(
            np.linspace(
                self.input["temperature_start"],
                self.input["temperature_end"],
                int(self.input["temperature_steps"])
            ),
            np.array(volume_lst)[arg_lst],
        )

        self._output["volumes"] = volume_mesh
        self._output["temperatures"] = temperature_mesh
        with self.project_hdf5.open("output") as hdf5_out:
            for key, val in self._output.items():
                hdf5_out[key] = val

    def optimise_volume(self, bulk_eng):
        v0_lst, free_eng_lst, entropy_lst, cv_lst = [], [], [], []
        for i, [t, free_energy, cv, entropy, v] in enumerate(
                zip(self["output/temperatures"].T,
                    self["output/free_energy"].T,
                    self["output/cv"].T,
                    self["output/entropy"].T,
                    self["output/volumes"].T)):
            fit = np.poly1d(np.polyfit(v, free_energy + bulk_eng, int(self.input["polynomial_degree"])))
            v0, ind = calc_v0_from_fit_funct(fit_funct=fit, x=v, save_range=0.0, return_ind=True)

            v0_lst.append(v0)
            free_eng_lst.append(fit([v0]))
            entropy_lst.append(entropy[ind])
            cv_lst.append(cv[ind])
        return v0_lst, free_eng_lst, entropy_lst, cv_lst
