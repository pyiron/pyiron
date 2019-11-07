# coding: utf-8
import numpy as np
from pyiron.base.master.parallel import JobGenerator


class USJobGenerator(JobGenerator):
    @property
    def parameter_list(self):
        """

        Returns:
            (list)
        """
        # For now no different kappa for different locs implementation!
        
        assert isinstance(self._job.input['grid_range'], list) or isinstance(self._job.input['grid_range'], np.ndarray)
        for loc in self._job.input['grid_range']:
            structure = self.search_structure(job.input['md_job'],job.input['cv_f'],loc)
            parameter_lst.append([np.round(loc,5), structure])
        return parameter_lst
    
    @staticmethod
    def search_structure(job,f,cv):
        # IMPLEMENT 

    @staticmethod
    def job_name(parameter):
        return "strain_" + str(parameter[0]).replace('.', '_')

    def modify_job(self, job, parameter):
        job.structure = parameter[1]
        return job
    
    
class US(AtomisticParallelMaster):
    def __init__(self, project, job_name='us'):
        """

        Args:
            project:
            job_name:
        """
        super(US, self).__init__(project, job_name)
        self.__name__ = 'us'
        self.__version__ = '0.1.0'

        # define default input
        self.input['grid_range'] = (np.linspace(0,1,10), 'cv grid')
        self.input['md_job'] = (None, 'job with md data for structure generation at cv grid points')
        self.input['cv_f'] = (None, 'function to calculate CV from job object')
        self._job_generator = USJobGenerator(self)

    def list_structures(self):
        if self.ref_job.structure is not None:
            return [parameter[1] for parameter in self._job_generator.parameter_list]
        else:
            return []

    # Fix these functions    
        
    def collect_output(self):
        if self.server.run_mode.interactive:
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
                print('job_id: ', job_id, ham.status)
                energy = ham["output/generic/energy_tot"][-1]
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
        if self.input['fit_type'] == "polynomial":
            self.fit_polynomial(fit_order=self.input['fit_order'])
        else:
            self._fit_eos_general(fittype=self.input['fit_type'])

    def get_structure(self, iteration_step=-1):
        """

        Returns: Structure

        """
        if not (self.structure is not None):
            raise AssertionError()
        if iteration_step == -1:
            snapshot = self.structure.copy()
            old_vol = snapshot.get_volume()
            new_vol = self["output/equilibrium_volume"]
            k = (new_vol / old_vol) ** (1. / 3.)
            new_cell = snapshot.cell * k
            snapshot.set_cell(new_cell, scale_atoms=True)
            return snapshot
        elif iteration_step == 0:
            return self.structure
        else:
            raise ValueError('iteration_step should be either 0 or -1.')
