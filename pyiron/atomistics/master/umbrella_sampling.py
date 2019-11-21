# coding: utf-8
import numpy as np
import matplotlib.pyplot as pt
from molmod.units import *
import subprocess, os

from pyiron.atomistics.master.parallel import AtomisticParallelMaster
from pyiron.base.master.parallel import JobGenerator
from pyiron.atomistics.structure.atoms import Atoms

class USJobGenerator(JobGenerator):
    @property
    def parameter_list(self):
        '''

        Returns:
            (list)
        '''
        # Check if all post processing parameters are correctly defined
        assert self._job.input['h_min'] is not None
        assert self._job.input['h_max'] is not None
        assert self._job.input['h_bins'] is not None

        # Create parameter list
        parameter_lst = []
        for (loc,structure) in zip(self._job.input['cv_grid'],self._job.structures):
            parameter_lst.append([np.round(loc,5), structure])
        return parameter_lst

    @staticmethod
    def job_name(parameter):
        if isinstance(parameter[0], list) or isinstance(parameter[0], np.ndarray):
            return 'us_' + '__'.join([str(loc).replace('.', '_').replace('-', 'm') for loc in parameter[0]])
        else:
            return 'us_' + str(parameter[0]).replace('.', '_').replace('-', 'm')

    def modify_job(self, job, parameter):
        # For now, no different kappa for different locs implementation!
        job.input['temp'] = self._job.input['temp']
        job.structure = parameter[1]
        job.set_us(self._job.input['cvs'], self._job.input['kappa'], parameter[0], fn_colvar='COLVAR', stride=self._job.input['stride'], temp=self._job.input['temp'])
        return job


class US(AtomisticParallelMaster):
    def __init__(self, project, job_name='us'):
        '''

        Args:
            project:
            job_name:
        '''
        super(US, self).__init__(project, job_name)
        self.__name__ = 'us'
        self.__version__ = '0.1.0'

        # define default input
        self.input['kappa']       = (1.*kjmol, 'force constant of the harmonic bias potential')
        self.input['stride']      = (10, 'step for output printed to COLVAR output file.')
        self.input['temp']        = (300*kelvin, 'the system temperature')
        self.input['cv_grid']     = (list(np.linspace(0,1,10)), 'cv grid, has to be a list')
        self.input['cvs']         = ([('distance', [0,1])], 'cv(s), see Yaff input for description')

        self.input['h_min']       = (None , 'lowest value(s) of the cv(s) for WHAM')
        self.input['h_max']       = (None , 'highest value(s) of the cv(s) for WHAM')
        self.input['h_bins']      = (None , 'bins between h_min and h_max for WHAM')

        self.input['periodicity'] = (None , 'periodicity of cv(s)')
        self.input['tol']         = (0.00001 , 'WHAM converges if free energy changes < tol')

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

        cv = cv_f(job).reshape(-1,len(self.input['cvs']))
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


    def check_overlap(self):
        '''
            Checks overlap between different umbrella simulations. Only works for 1D!

            **Arguments**
        '''
        pt.figure()
        for job_id in self.child_ids:
            job = self.project_hdf5.inspect(job_id)
            pt.plot(job['output/enhanced/cv'])
        pt.show()


    def wham(self, h_min, h_max, bins, f_metadata, f_fes, periodicity=None, tol=0.00001):
        '''
            Performs the weighted histogram analysis method to calculate the free energy surface

            **Arguments**

            h_min   lowest value that is taken into account, float or list if more than one cv is biased

            h_max   highest value that is taken into account, float or list if more than one cv is biased
                    if one whole trajectory is outside of these borders an error occurs

            bins    number of bins between h_min and h_max

            periodicity
                    periodicity of the collective variable
                    1D: either a number, 'pi' for angles (2pi periodic) or an empty string ('') for periodicity of 360 degrees
                    2D: either a number, 'pi' for angles (2pi periodic) or 0 if no periodicity is required

            tol     if no free energy value changes between iteration for more than tol, wham is converged
        '''


        if isinstance(h_min,(int,float)):
            cmd = os.path.join(self.project.root_path,'wham')
            cmd += ' '
            if not periodicity is None:
                cmd += 'P{} '.format(periodicity)
            cmd += ' '.join(map(str,[h_min,h_max,int(bins),tol,self.input['temp'],0,f_metadata,f_fes]))

        elif isinstance(h_min,list) and len(h_min) == 2:
            cmd =  os.path.join(self.project.root_path,'wham-2d')
            cmd += ' '
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


    def get_structure(self, iteration_step=-1):
        '''

        Returns: Structure at free energy minimum

        '''

        # Read minimal energy from fes
        # Read corresponding job
        # return average structure

        raise NotImplementedError()

    def to_hdf(self, hdf=None, group_name=None):
        super(US, self).to_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open('input') as hdf5_input:
            self.input.to_hdf(hdf5_input)

        with self.project_hdf5.open('input/structures') as hdf5_input:
            for n,(loc,structure) in enumerate(zip(self.input['cv_grid'],self.structures)):
                #name = str(loc).replace('.', '_').replace('-', 'm') if isinstance(loc,(int,float)) else ','.join([str(l).replace('.', '_').replace('-', 'm') for l in loc])
                name = 'cv' + str(n)
                self.structure.to_hdf(hdf5_input,group_name=name)

    def from_hdf(self, hdf=None, group_name=None):
        super(US, self).from_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open('input') as hdf5_input:
            self.input.from_hdf(hdf5_input)

        self.structures = []
        with self.project_hdf5.open('input/structures') as hdf5_input:
            for n,loc in enumerate(self.input['cv_grid']):
                #name = str(loc).replace('.', '_').replace('-', 'm') if isinstance(loc,(int,float)) else ','.join([str(l).replace('.', '_').replace('-', 'm') for l in loc])
                name = 'cv' + str(n)
                self.structures.append(Atoms().from_hdf(hdf5_input,group_name=name))

    def collect_output(self):
        def convert_val(val,unit=None):
            scale = 1 if unit is None else unit
            if isinstance(val, list) or isinstance(val, np.ndarray):
                return [str(l/scale) for l in val]
            else:
                return str(val/scale)

        # Files to store data of wham
        f_metadata = os.path.join(self.working_directory, 'metadata')
        f_fes      = os.path.join(self.working_directory, 'fes.dat')


        with open(f_metadata, 'w') as f:
            for job_id in self.child_ids:
                job = self.project_hdf5.inspect(job_id)
                print('job_id: ', job_id, job.status)
                loc = convert_val(job['input/generic/enhanced/loc'])
                kappa = convert_val(job['input/generic/enhanced/kappa'],unit=kjmol)
                f.write('{}/COLVAR\t'.format(job.working_directory) + '\t'.join(loc) + '\t' + '\t'.join(kappa) + '\n') # format of colvar needs to be TIME CV1 (CV2)

        # Execute wham code
        self.wham(self.input['h_min'], self.input['h_max'], self.input['h_bins'], f_metadata, f_fes, periodicity=self.input['periodicity'], tol=self.input['tol'])

        # Process output of wham code
        data = np.loadtxt(os.path.join(self.working_directory,'fes.dat'))

        if len(loc) == 1:
            bins = data[:,0]
            fes = data[:,1]
        elif len(loc) == 2:
            bins = data[:,0:2]
            fes = data[:,2]

        with self.project_hdf5.open('output') as hdf5_out:
            hdf5_out['bins'] = bins
            hdf5_out['fes'] = fes
