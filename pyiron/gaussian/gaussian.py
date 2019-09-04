# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import os,subprocess
import numpy as np

from pyiron.base.generic.parameters import GenericParameters
from pyiron.atomistics.structure.atoms import Atoms
from pyiron.atomistics.job.atomistic import AtomisticGenericJob,GenericOutput,Trajectory
from molmod.io.fchk import FCHKFile
from molmod.units import amu,angstrom,electronvolt
from molmod.periodic import periodic
import tamkin

__author__ = "Jan Janssen, Sander Borgmans"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - " \
                "- Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = ""
__email__ = ""
__status__ = "trial"
__date__ = "Aug 27, 2019"


class Gaussian(AtomisticGenericJob):
    def __init__(self, project, job_name):
        super(Gaussian, self).__init__(project, job_name)
        self.__name__ = "Gaussian"
        self._executable_activate(enforce=True)
        self.input = GaussianInput()
        self.output = GaussianOutput(job=self)

    def write_input(self):
        input_dict = {'mem': self.server.memory_limit,
                      'verbosity': self.input['verbosity'],
                      'lot': self.input['lot'],
                      'basis_set': self.input['basis_set'],
                      'jobtype' : self.input['jobtype'],
                      'settings' : self.input['settings'],
                      'title' : self.input['title'],
                      'spin_mult': self.input['spin_mult'],
                      'charge': self.input['charge'],
                      'symbols': self.structure.get_chemical_symbols().tolist(),
                      'pos': self.structure.positions}
        write_input(input_dict=input_dict, working_directory=self.working_directory)

    def collect_output(self):
        output_dict = collect_output(output_file=os.path.join(self.working_directory, 'input.fchk'))
        with self.project_hdf5.open("output") as hdf5_output:
            for k, v in output_dict.items():
                hdf5_output[k] = v

    def to_hdf(self, hdf=None, group_name=None):
        super(Gaussian, self).to_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf5_input:
            self.structure.to_hdf(hdf5_input)
            self.input.to_hdf(hdf5_input)

    def from_hdf(self, hdf=None, group_name=None):
        super(Gaussian, self).from_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf5_input:
            self.input.from_hdf(hdf5_input)
            self.structure = Atoms().from_hdf(hdf5_input)

    def visualize_MO(self,index):
        n_MO = self.output.scf_density.shape[0]
        assert index > 0 and index <= n_MO
        assert len(self.output.numbers) < 50 # check whether structure does not become too large for interactive calculation of cube file

        #make cube file
        path = self.path+'_hdf5/'+self.name+'/input'
        out = subprocess.check_output(
                "ml load Gaussian/g16_E.01-intel-2018a; cubegen 1 MO={} {}.fchk {}.cube".format(index,path,path),
                stderr=subprocess.STDOUT,
                universal_newlines=True,
                shell=True,
            )
        #visualize cube file
        try:
            import nglview
        except ImportError:
            raise ImportError("The animate_nma_mode() function requires the package nglview to be installed")

        picture = nglview.show_ase(self.structure)
        picture.add_component('{}.cube'.format(path), opacity=.3)
        return picture

    def do_nma(self):
        mol = tamkin.Molecule(self.output.numbers, self.output.positions, self.output.masses, self.output.energy_tot, self.output.forces *-1, self.output.hessian)
        self.nma = tamkin.NMA(mol)

    def animate_nma_mode(self,index,amplitude=1.0,frames=24,spacefill=False,particle_size=0.5):
        coordinates = self.nma.coordinates
        symbols = [periodic[n].symbol for n in self.nma.numbers]

        mode = self.nma.modes[:,index]
        if self.nma.masses3 is not None:
            mode /= np.sqrt(self.nma.masses3)
        mode /= np.linalg.norm(mode)

        positions = np.zeros((frames,len(symbols),3))

        for frame in range(frames):
            factor = amplitude*np.sin(2*np.pi*float(frame)/frames)
            positions[frame] = coordinates + factor*mode.reshape((-1,3))

        try:
            import nglview
        except ImportError:
            raise ImportError("The animate_nma_mode() function requires the package nglview to be installed")

        animation = nglview.show_asetraj(Trajectory(positions,self.structure))
        if spacefill:
            animation.add_spacefill(radius_type='vdw', scale=0.5, radius=particle_size)
            animation.remove_ball_and_stick()
        else:
            animation.add_ball_and_stick()
        return animation

class GaussianInput(GenericParameters):
    def __init__(self, input_file_name=None):
        super(GaussianInput, self).__init__(input_file_name=input_file_name, table_name="input_inp", comment_char="#")

    def load_default(self):
        """
        Loading the default settings for the input file.
        """
        input_str = """\
lot HF
basis_set 6-311G(d,p)
spin_mult 1
charge 0
"""
        self.load_string(input_str)

class GaussianOutput(GenericOutput):
    """
    Handles the output from a Gaussian simulation.
    Adds extra properties to handle for NMA
    """

    @property
    def numbers(self):
        return self._job['output/generic/numbers']

    @property
    def masses(self):
        return self._job['output/generic/masses']

    @property
    def hessian(self):
        return self._job['output/generic/hessian']

    @property
    def scf_density(self):
        return self._job['output/generic/dft/scf_density']


def write_input(input_dict,working_directory='.'):
    # Comments can be written with ! in Gaussian
    # Load dictionary
    lot          = input_dict['lot']
    basis_set    = input_dict['basis_set']
    spin_mult    = input_dict['spin_mult'] # 2S+1
    charge       = input_dict['charge']
    symbols      = input_dict['symbols']
    pos          = input_dict['pos']
    assert pos.shape[0] == len(symbols)

    # Optional elements
    if not input_dict['mem'] is None:
        mem = input_dict['mem'] + 'B' * (input_dict['mem'][-1]!='B') # check if string ends in bytes
    else:
        mem = "800MB" # default allocation

    if not input_dict['jobtype'] is None:
        jobtype = input_dict['jobtype']
    else:
        jobtype = "" # corresponds to sp

    if not input_dict['title'] is None:
        title = input_dict['title']
    else:
        title = "no title"

    if not input_dict['settings'] is None:
        settings = input_dict['settings'] # dictionary {key: [options]}
    else:
        settings = {}

    verbosity_dict={'low':'t','normal':'n','high':'p'}
    if not input_dict['verbosity'] is None:
        verbosity  = input_dict['verbosity']
        if verbosity in verbosity_dict:
            verbosity = verbosity_dict[verbosity]
    else:
        verbosity='n'

    # Parse settings
    settings_string = ""
    for key,valuelst in settings.items():
        if not isinstance(valuelst, list):
            valuelst = [valuelst]
        option = key + "({}) ".format(",".join(valuelst))*(len(valuelst)>0)
        settings_string += option

    # Write to file
    route_section = "#{} {}/{} {} {}\n\n".format(verbosity,lot,basis_set,jobtype,settings_string)
    with open(os.path.join(working_directory, 'input.com'), 'w') as f:
        f.write("%mem={}\n".format(mem))
        f.write("%chk=input.chk\n")
        f.write(route_section)
        f.write("{}\n\n".format(title))
        f.write("{} {}\n".format(charge,spin_mult))
        for n,p in enumerate(pos):
            f.write(" {}\t{: 1.6f}\t{: 1.6f}\t{: 1.6f}\n".format(symbols[n],p[0],p[1],p[2]))
        f.write('\n\n') # don't know whether this is still necessary in G16


# we could use theochem iodata, should be more robust than molmod.io
# but we require the latest iodata for this, not the conda version
def fchk2dict(fchk):
    # probably still some data missing
    # check job type, for now implement basics (SP=single point, FOpt = full opt, Freq = frequency calculation)
    if not fchk.command.lower() in ['sp','fopt','freq']:
        raise NotImplementedError

    # Basic information
    fchkdict = {}
    fchkdict['jobtype']     = fchk.command.lower()
    fchkdict['lot']         = fchk.lot
    fchkdict['basis_set']   = fchk.basis

    fchkdict['generic/numbers']     = fchk.fields.get('Atomic numbers')
    fchkdict['generic/masses']      = fchk.fields.get("Real atomic weights")*amu
    fchkdict['generic/charges']     = fchk.fields.get('Mulliken Charges')
    fchkdict['generic/dipole']      = fchk.fields.get('Dipole Moment')
    fchkdict['generic/dft/n_electrons']         = fchk.fields.get('Number of electrons')
    fchkdict['generic/dft/n_alpha_electrons']   = fchk.fields.get('Number of alpha electrons')
    fchkdict['generic/dft/n_beta_electrons']    = fchk.fields.get('Number of beta electrons')
    fchkdict['generic/dft/n_basis_functions']   = fchk.fields.get('Number of basis functions')

    # Orbital information
    fchkdict['generic/dft/alpha_orbital_e']     = fchk.fields.get('Alpha Orbital Energies')
    fchkdict['generic/dft/beta_orbital_e']      = fchk.fields.get('Beta Orbital Energies')

    # Densities
    fchkdict['generic/dft/scf_density']         = _triangle_to_dense(fchk.fields.get('Total SCF Density'))
    fchkdict['generic/dft/spin_scf_density']    = _triangle_to_dense(fchk.fields.get('Spin SCF Density'))

    if fchk.lot.upper() in ['MP2', 'MP3', 'CC', 'CI']:
        # only one of the lots should be present, hence using the same key
        fchkdict['generic/dft/post_scf_density']      = _triangle_to_dense(fchk.fields.get('Total {} Density'.format(fchk.lot)))
        fchkdict['generic/dft/post_spin_scf_density'] = _triangle_to_dense(fchk.fields.get('Spin {} Density'.format(fchk.lot)))

    # Specific job information
    if fchkdict['jobtype'] == 'fopt':
        fchkdict['generic/positions']     = fchk.get_optimization_coordinates()/angstrom
        fchkdict['generic/energy_tot']    = fchk.get_optimization_energies()/electronvolt
        fchkdict['generic/forces']        = fchk.get_optimization_gradients()/(electronvolt/angstrom) * -1

    if fchkdict['jobtype'] == 'freq':
        fchkdict['generic/positions']     = fchk.fields.get('Current cartesian coordinates').reshape([-1, 3])/angstrom
        fchkdict['generic/forces']        = fchk.fields.get('Cartesian Gradient').reshape([-1, 3])/(electronvolt/angstrom) *-1
        fchkdict['generic/hessian']       = fchk.get_hessian()/(electronvolt**2/angstrom**2)
        fchkdict['generic/energy_tot']    = fchk.fields.get('Total Energy')/electronvolt

    if fchkdict['jobtype'] == 'sp':
        fchkdict['generic/positions']     = fchk.fields.get('Current cartesian coordinates').reshape([-1, 3])/angstrom
        fchkdict['generic/energy_tot']    = fchk.fields.get('Total Energy')/electronvolt

    return fchkdict

def collect_output(output_file):
    '''
    # First check log file if terminated ok
    with open(log_file,'r'):
        lines = f.readlines()
    if not 'Normal termination of Gaussian' in lines[-1]: raise ValueError('No normal termination of Gaussian')
    '''
    # Read output
    fchk = FCHKFile(output_file)

    # Translate to dict
    output_dict = fchk2dict(fchk)

    return output_dict

# function from theochem iodata
def _triangle_to_dense(triangle):
    """Convert a symmetric matrix in triangular storage to a dense square matrix.
    Parameters
    ----------
    triangle
        A row vector containing all the unique matrix elements of symmetric
        matrix. (Either the lower-triangular part in row major-order or the
        upper-triangular part in column-major order.)
    Returns
    -------
    ndarray
        a square symmetric matrix.
    """
    if triangle is None: return None
    nrow = int(np.round((np.sqrt(1 + 8 * len(triangle)) - 1) / 2))
    result = np.zeros((nrow, nrow))
    begin = 0
    for irow in range(nrow):
        end = begin + irow + 1
        result[irow, :irow + 1] = triangle[begin:end]
        result[:irow + 1, irow] = triangle[begin:end]
        begin = end
    return result
