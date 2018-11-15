import numpy as np
from pyiron.atomistics.hamilton.murnaghan import Murnaghan


class MurnaghanInt(Murnaghan):
    def create_jobs(self):
        if self.ref_job.server.run_mode.interactive:
            ham = self.ref_job.copy()
            basis_ref = ham.structure

            vol_min = 1 - self.input['vol_range']
            vol_max = 1 + self.input['vol_range']

            self.submission_status.submitted_jobs = 0
            ham.master_id = self.job_id
            for strain in np.linspace(vol_min, vol_max, self.input['num_points']):
                strain = np.round(strain, 7)  # prevent too long job_names
                self._logger.debug("Murnaghan child project {}".format(self.project.__str__()))
                basis = basis_ref.copy()
                basis.set_cell(basis.cell * strain ** (1. / 3.), scale_atoms=True)
                ham.structure = basis
                ham.run()
                self._logger.info('Murnaghan: finished job {}'.format(ham.job_name))
                self.ref_job.structure = basis_ref
            self.structure = basis_ref
            ham.interactive_close()
            self.status.collect = True
            self.run()
        else:
            super(MurnaghanInt, self).create_jobs()

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