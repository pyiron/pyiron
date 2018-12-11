# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
from pyiron.atomistics.master.murnaghan import Murnaghan

__author__ = "Jan Janssen"
__copyright__ = "Copyright 2017, Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


class MurnaghanInt(Murnaghan):
    def run_static(self):
        if self.ref_job.server.run_mode.interactive:
            ham = self.ref_job.copy()
            ham.master_id = self.job_id

            for strain in self._job_generator.parameter_list:
                ham = self._job_generator.modify_job(job=ham, parameter=strain)
                ham.run()

            ham.interactive_close()
            self.status.collect = True
            self.run()
        else:
            super(MurnaghanInt, self).run_static()

    def collect_output(self):
        erg_lst, vol_lst, err_lst, id_lst = [], [], [], []
        for job_id in self.child_ids:
            ham = self.project_hdf5.inspect(job_id)
            print('job_id: ', job_id, ham.status)
            energy = ham["output/generic/energy_tot"]
            volume = ham["output/generic/volume"]
            erg_lst = [np.mean(eng) for eng in energy]
            err_lst = [np.var(eng) for eng in energy]
            vol_lst = volume
            id_lst.append(job_id)
        vol_lst = np.array(vol_lst)
        erg_lst = np.array(erg_lst)
        err_lst = np.array(err_lst)
        arg_lst = np.argsort(vol_lst)

        self._output["volume"] = vol_lst[arg_lst]
        self._output["energy"] = erg_lst[arg_lst]
        self._output["error"] = err_lst[arg_lst]