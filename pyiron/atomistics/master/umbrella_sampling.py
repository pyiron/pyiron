# coding: utf-8
import numpy as np
from molmod.units import *
import subprocess, os

from pyiron.atomistics.master.parallel import AtomisticParallelMaster
from pyiron.base.master.parallel import JobGenerator



class USJobGenerator(JobGenerator):
    @property
    def parameter_list(self):
        """

        Returns:
            (list)
        """
        # For now no different kappa for different locs implementation!
        parameter_lst = []
        assert isinstance(self._job.input['cv_grid'], list) or isinstance(self._job.input['cv_grid'], np.ndarray)
        for (loc,structure) in zip(self._job.input['cv_grid'],self._job.structures):
            parameter_lst.append([np.round(loc,3), structure])
        return parameter_lst

    @staticmethod
    def job_name(parameter):
        if isinstance(parameter[0], list) or isinstance(parameter[0], np.ndarray):
            return 'us_' + '__'.join([str(loc).replace('.', '_').replace('-', 'm') for loc in parameter[0]])
        else:
            return 'us_' + str(parameter[0]).replace('.', '_').replace('-', 'm')

    def modify_job(self, job, parameter):
        job.input['temp'] = self._job.input['temp']
        job.structure = parameter[1]
        job.set_us(self._job.input['ics'], self._job.input['kappa'], parameter[0], fn_colvar='COLVAR', stride=self._job.input['stride'], temp=self._job.input['temp'])
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
        self.input['kappa']      = (1.*kjmol, 'the value of the force constant of the harmonic bias potential')
        self.input['stride']     = (10, 'the number of steps after which the internal coordinate values and bias are printed to the COLVAR output file.')
        self.input['temp']       = (300*kelvin, 'the system temperature')

        self.input['cv_grid']    = (list(np.linspace(0,1,10)), 'cv grid, has to be a list')
        self.input['ics']        = ([('distance', [0,1])], 'ics')

        self.structures = None   # list with structures corresponding to grid points
        self._job_generator = USJobGenerator(self)

    def list_structures(self):
        return self.structures

    def generate_structures_traj(self,job,cv_f):
        '''
            Generates structure list based on cv grid and cv function using the trajectory data from another job (e.g. MD or MTD job)

            **Arguments**

            job      job object which contains enough snapshots in the region of interest
            cv_f     function object that takes a job object as input and returns the corresponding CV(s) list
        '''

        cv = cv_f(job).reshape(-1,len(self.input['ics']))
        idx = np.zeros(len(self.input['cv_grid']),dtype=int)
        for n,loc in enumerate(self.input['cv_grid']):
            idx[n] = np.argmin(np.linalg.norm(loc-cv,axis=-1))

        return [job.get_structure(i) for i in idx]

    def generate_structures_ref(self,f):
        '''
            Generates structure list based on cv grid and reference structure

            **Arguments**

            f     function object that takes the reference structure object and a cv change as input and returns the altered structure
        '''

        assert self.ref_job.structure is not None
        structures = []
        for loc in self.input['cv_grid']:
            structures.append(f(self.ref_job.structure,loc))
        return structures

    def collect_output(self):
        if self.server.run_mode.interactive:
            job = self.project_hdf5.inspect(self.child_ids[0])
            data = np.loadtxt(os.path.join(job.working_directory,'COLVAR'))
            time_lst = data[:,0]
            cv_lst = data[:,1]
            bias_lst = data[:,2]
            self._output["time"] = time_lst
            self._output["cv"] = cv_lst
            self._output["bias"] = bias_lst
        else:
            time_lst, cv_lst, bias_lst, id_lst = [], [], [], []
            for job_id in self.child_ids:
                job = self.project_hdf5.inspect(job_id)
                print('job_id: ', job_id, job.status)
                data = np.loadtxt(os.path.join(job.working_directory,'COLVAR'))
                time = data[:,0]
                cv = data[:,1]
                bias = data[:,2]
                time_lst.append(time)
                cv_lst.append(cv)
                bias_lst.append(bias)
                id_lst.append(job_id)
            time_lst = np.array(time_lst)
            cv_lst = np.array(cv_lst)
            bias_lst = np.array(bias_lst)
            id_lst = np.array(id_lst)

            self._output["time"] = time_lst
            self._output["cv"] = cv_lst
            self._output["bias"] = bias_lst
            self._output["id"] = id_lst

        with self.project_hdf5.open("output") as hdf5_out:
            for key, val in self._output.items():
                hdf5_out[key] = val

    def check_overlap(self):
        '''
            Checks overlap between different umbrella simulations. Only works for 1D!

            **Arguments**
        '''
        pt.figure()
        for job_id in self.child_ids:
            pt.plot(job['output/cv'])
        pt.show()

    def wham(self, h_min, h_max, bins, periodicity=None, tol=0.00001):
        '''
            Performs the weighted histogram analysis method to calculate the free energy surface

            **Arguments**

            h_min   lowest value that is taken into account, float or list if more than one cv is biased

            h_max   highest value that is taken into account, float or list if more than one cv is biased
                    if one whole trajectory is outside of these borders an error occurs

            bins    number of bins between h_min and h_max, int

            periodicity
                    periodicity of the collective variable
                    1D: either a number, 'pi' for angles (2pi periodic) or an empty string ('') for periodicity of 360 degrees
                    2D: either a number, 'pi' for angles (2pi periodic) or 0 if no periodicity is required

            tol     if no free energy value changes between iteration for more than tol, wham is converged
        '''

        def convert_val(val,unit=None):
            scale = 1 if unit is None else unit
            if isinstance(val, list) or isinstance(val, np.ndarray):
                return [str(l/scale) for l in val]
            else:
                return str(val/scale)

        f_metadata = os.path.join(self.working_directory, 'metadata')
        f_fes      = os.path.join(self.working_directory, 'fes.dat')

        with open(f_metadata, 'w') as f:
            for job_id in self.child_ids:
                job = self.project_hdf5.inspect(job_id)
                print('job_id: ', job_id, job.status)
                loc = convert_val(job['input/generic/enhanced/loc'])
                kappa = convert_val(job['input/generic/enhanced/kappa'],unit=kjmol)
                f.write('{}/COLVAR\t'.format(job.working_directory) + '\t'.join(loc) + '\t' + '\t'.join(kappa) + '\n') # format of colvar needs to be TIME CV1 (CV2)

        if len(loc) == 1:
            cmd = './wham '
            if not periodicity is None:
                cmd += 'P{} '.format(periodicity)
            cmd += ' '.join(map(str,[h_min,h_max,int(bins),tol,self.input['temp'],0,f_metadata,f_fes]))

        elif len(loc) == 2:
            cmd = './wham-2d '
            periodic = ['Px='+str(periodicity[0]) if not periodicity[0] is None else '0', 'Py='+str(periodicity[1]) if not periodicity[1] is None else '0']
            for i in range(2):
                cmd += ' '.join([periodic[i],h_min[i],h_max[i],int(bins[i])])
            cmd += ' '.join([tol,self.input['temp'],0,f_metadata,f_fes,1])
        else:
            raise NotImplementedError()

        subprocess.check_output(
            cmd,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            shell=True
        )

        data = np.loadtxt(os.path.join(self.working_directory,'fes.dat'))

        if len(loc) == 1:
            bins = data[:,0]
            fes = data[:,1]
        elif len(loc) == 2:
            bins = data[:,0:2]
            fes = data[:,2]

        with self.project_hdf5.open("output") as hdf5_out:
            hdf5_out["bins"] = bins
            hdf5_out["fes"] = fes

    def get_structure(self, iteration_step=-1):
        """

        Returns: Structure at free energy minimum

        """

        # Read minimal energy from fes
        # Read corresponding job
        # return average structure

        raise NotImplementedError()
