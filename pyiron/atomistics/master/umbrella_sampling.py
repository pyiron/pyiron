# coding: utf-8
import numpy as np
from molmod.units import *
import subprocess

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

    """
    def collect_output(self):
        '''
            Executes the plumed post-processing functions

        '''
        
        ickinds   = np.array([ic[0] for ic in self.input['ics']],dtype='S22')
        icindices = np.array([np.array(ic[1])+1 for ic in self.input['ics']]) # plumed starts counting from 1
        
        locs = []
        for job_id in self.child_ids:
            job = self.project_hdf5.inspect(job_id)
            print('job_id: ', job_id, job.status)
            loc = job['input/generic/enhanced/loc']
            if isinstance(loc, list) or isinstance(loc, np.ndarray):
                locs.append(",".join([str(l) for l in loc]))
            else:
                locs.append(str(loc))
            us.load(job_id).write_traj(os.path.join(self.working_directory, 'alltraj.xyz'), append=True)

        with open(os.path.join(self.working_directory, 'plumed.dat'), 'w') as f:
            #set units to atomic units
            f.write('UNITS LENGTH=Bohr ENERGY=kj/mol TIME=atomic \n')
            #define ics
            for i, kind in enumerate(ickinds):
                if isinstance(kind, bytes):
                    kind = kind.decode()
                if len(icindices[i] > 0):
                    f.write('ic%i: %s ATOMS=%s \n' %(i, kind.upper(), ','.join([str(icidx) for icidx in icindices[i]])))
                else:
                    f.write('ic%i: %s \n' %(i, kind.upper(), ','.join([str(icidx) for icidx in icindices[i]])))
                f.write('umbrella: RESTRAINT ARG=%s KAPPA=%s AT=@replicas:{\n %s \n} \n' %(
                ','.join([ 'ic%i' %i for i in range(len(ickinds))]),
                kappa/kjmol, '\n'.join([str(loc) for loc in locs])
            ))

            # Current implementation only works for 1D umbrella sampling
            f.write('hh: WHAM_HISTOGRAM ARG=%s BIAS=umbrella.bias TEMP=%s GRID_MIN=%s GRID_MAX=%s GRID_BIN=%s \n' %(
                ','.join([ 'ic%i' %i for i in range(len(ickinds))]), self.input['temp'],
                np.min(locs), np.max(locs), len(locs)
            ))

            f.write('fes: CONVERT_TO_FES GRID=hh TEMP=%s \n' %(self.input['temp']))
            f.write('DUMPGRID GRID=fes FILE=fes.dat ')

            subprocess.check_output(
            'ml load PLUMED/2.5.2-intel-2019a-Python-3.7.2; mpirun -np 6 plumed driver --ixyz alltraj.xyz --multi 6',
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            shell=True
        )
    """
        
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
        def convert_val(val):
            if isinstance(val, list) or isinstance(val, np.ndarray):
                return [str(l) for l in val]
            else:
                return str(val)
        
        f_metadata = os.path.join(self.working_directory, 'metadata')
        f_fes      = os.path.join(self.working_directory, 'fes.dat')
        
        with open(f_metadata, 'w') as f:
            for job_id in self.child_ids:
                job = self.project_hdf5.inspect(job_id)
                print('job_id: ', job_id, job.status)
                loc = convert_val(job['input/generic/enhanced/loc'])
                kappa = convert_val(job['input/generic/enhanced/kappa'])
                f.write('{}/COLVAR\t'.format(job.working_directory) + '\t'.join(loc) + '\t' + '\t'.join(loc) + '\n') # format of colvar needs to be TIME CV1 (CV2)
        
        if len(loc) == 1:
            cmd = 'wham '
            if not periodicity is None:
                cmd += 'P{} '.format(periodicity)
            cmd += ' '.join([hmin,hmax,int(bins),tol,self.input['temp'],0,f_metadata,f_fes])
            
        elif len(loc) == 2:
            cmd = 'wham-2d '
            periodic = ['Px='+str(periodicity[0]) if not periodicity[0] is None else '0', 'Py='+str(periodicity[1]) if not periodicity[1] is None else '0']
            for i in range(2):
                cmd += ' '.join([periodic[i],hmin[i],hmax[i],int(bins[i])])
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
        """

        Returns: Structure at free energy minimum

        """

        # Read minimal energy from fes
        # Read corresponding job
        # return average structure

        raise NotImplementedError()
