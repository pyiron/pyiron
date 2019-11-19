# coding: utf-8
import numpy as np
from molmod.units import *
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

        assert isinstance(self._job.input['cv_grid'], list) or isinstance(self._job.input['cv_grid'], np.ndarray)
        for (loc,struc) in zip(self._job.input['cv_grid'],self._job.input['structures']):
            parameter_lst.append([np.round(loc,5), structure])
        return parameter_lst

    @staticmethod
    def job_name(parameter):
        if isinstance(parameter[0], list) or isinstance(parameter[0], np.ndarray):
            return 'us_' + '__'.join([str(loc).replace('.', '_') for loc in parameter[0]])
        else:
            return 'us_' + str(parameter[0]).replace('.', '_')

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

        self.input['cv_grid']    = (np.linspace(0,1,10), 'cv grid')
        self.input['structures'] = (None, 'list with structures corresponding to grid points')
        self.input['ics']        = ([('distance', [0,1])], 'ics')

        self._job_generator = USJobGenerator(self)

    def list_structures(self):
        return self.input['structures']

    def generate_structures_traj(self,job,cv_f):
        '''
            Generates structure list based on cv grid and cv function using the trajectory data from another job (e.g. MD or MTD job)

            **Arguments**

            job      job object which contains enough snapshots in the region of interest
            cv_f     function object that takes a structure object as input and returns the corresponding CV(s)
        '''
        frames = job['output/generic/positions'].shape[0]
        cv = np.zeros((frames,len(self.input['ics'])))

        for i in range(frames):
            structure = job.get_structure(i)
            cv[i] = cv_f(structure)

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


    def wham(self):
        '''
            Executes the plumed post-processing functions

        '''
        
        ickinds   = np.array([ic[0] for ic in self.input['ics']],dtype='S22')
        icindices = np.array([np.array(ic[1])+1 for ic in self.input['ics']]) # plumed starts counting from 1
        
        locs = []
        for job_id in self.child_ids:
            job = self.project_hdf5.inspect(job_id)
            print('job_id: ', job_id, job.status)
            loc = job.enhanced['loc']
            if isinstance(loc, list) or isinstance(loc, np.ndarray):
                locs.append(",".join([str(l) for l in loc]))
            else:
                locs.append(str(loc))
            job.write_traj(os.path.join(self.working_directory, 'alltraj.xyz'), append=True)

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
                kappa, '\n'.join([str(loc) for loc in locs])
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


    def get_structure(self, iteration_step=-1):
        """

        Returns: Structure at free energy minimum

        """

        # Read minimal energy from fes
        # Read corresponding job
        # return average structure

        raise NotImplementedError()
