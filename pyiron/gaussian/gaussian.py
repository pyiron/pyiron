# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import os
import subprocess
import re
import pandas
import numpy as np
import matplotlib.pyplot as pt

from pyiron.dft.job.generic import GenericDFTJob
from pyiron_base import GenericParameters, ImportAlarm
from pyiron.atomistics.structure.atoms import Atoms

try:
    from molmod.io.fchk import FCHKFile
    from molmod.units import amu, angstrom, electronvolt, centimeter, kcalmol
    from molmod.constants import lightspeed
    from molmod.periodic import periodic
    import tamkin
    import_alarm = ImportAlarm()
except ImportError:
    import_alarm = ImportAlarm(
        "Gaussian relies on the molmod and tamkin packages, but these are unavailable. Please ensure your python "
        "environment contains these."
    )


__author__ = "Jan Janssen, Sander Borgmans"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - " \
                "- Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = ""
__email__ = ""
__status__ = "trial"
__date__ = "Aug 27, 2019"


class Gaussian(GenericDFTJob):
    @import_alarm
    def __init__(self, project, job_name):
        super(Gaussian, self).__init__(project, job_name)
        self.__name__ = "Gaussian"
        self._executable_activate(enforce=True)
        self.input = GaussianInput()

    def write_input(self):
        input_dict = {'mem': self.server.memory_limit,
                      'cores': self.server.cores,
                      'verbosity': self.input['verbosity'],
                      'lot': self.input['lot'],
                      'basis_set': self.input['basis_set'],
                      'jobtype' : self.input['jobtype'],
                      'settings' : self.input['settings'],
                      'title' : self.input['title'],
                      'spin_mult': self.input['spin_mult'],
                      'charge': self.input['charge'],
                      'bsse_idx': self.input['bsse_idx'],
                      'symbols': self.structure.get_chemical_symbols().tolist(),
                      'pos': self.structure.positions
                      }
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

    def log(self):
        with open(os.path.join(self.working_directory, 'input.log')) as f:
            print(f.read())

    def calc_minimize(self, electronic_steps=None, ionic_steps=None, algorithm=None, ionic_forces=None):
        """
            Function to setup the hamiltonian to perform ionic relaxations using DFT. The convergence goal can be set using
            either the iconic_energy as an limit for fluctuations in energy or the iconic_forces.

            **Arguments**

                algorithm: SCF algorithm
                electronic_steps (int): maximum number of electronic steps per electronic convergence
                ionic_steps (int): maximum number of ionic steps
                ionic_forces ('tight' or 'verytight'): convergence criterium for Berny opt (optional)
        """
        settings = {}
        opt_settings = []

        if electronic_steps is not None:
            if not 'SCF' in settings:
                settings['SCF'] = []
            settings['SCF'].append("MaxCycle={}".format(electronic_steps))

        if ionic_steps is not None:
            opt_settings.append("MaxCycles={}".format(ionic_steps))

        if algorithm is not None:
            if not 'SCF' in settings:
                settings['SCF'] = []
            settings['SCF'].append(algorithm)

        if ionic_forces is not None:
            assert isinstance(ionic_forces,str)
            opt_settings.append(ionic_forces)

        self.input['jobtype'] = 'opt' + '({})'.format(",".join(opt_settings))*(len(opt_settings)>0)
        if not isinstance(self.input['settings'],dict):
            self.input['settings'] = settings
        else:
            self.input['settings'].update(settings)

        super(Gaussian, self).calc_minimize(
            electronic_steps=electronic_steps,
            ionic_steps=ionic_steps,
            algorithm=algorithm,
            ionic_force_tolerance=ionic_forces
        )

    def calc_static(self, electronic_steps=None, algorithm=None):
        """
            Function to setup the hamiltonian to perform static SCF DFT runs

            **Arguments**

                algorithm (str): SCF algorithm
                electronic_steps (int): maximum number of electronic steps, which can be used to achieve convergence
        """
        settings = {}
        if electronic_steps is not None:
            if not 'SCF' in settings:
                settings['SCF'] = []
            settings['SCF'].append("MaxCycle={}".format(electronic_steps))

        if algorithm is not None:
            if not 'SCF' in settings:
                settings['SCF'] = []
            settings['SCF'].append(algorithm)

        self.input['jobtype'] = 'sp'
        if not isinstance(self.input['settings'],dict):
            self.input['settings'] = settings
        else:
            self.input['settings'].update(settings)

        super(Gaussian, self).calc_static(
            electronic_steps=electronic_steps,
            algorithm=algorithm
        )

    def calc_md(self, temperature=None,  n_ionic_steps=1000, time_step=None, n_print=100):
        raise NotImplementedError("calc_md() not implemented in Gaussian.")

    def print_MO(self):
        """
            Print a list of the MO's with the corresponding orbital energy and occupation.
        """

        n_MO = self.get('output/structure/dft/scf_density').shape[0]
        for n,index in enumerate(range(n_MO)):
            # print orbital information
            occ_alpha = int(self.get('output/structure/dft/n_alpha_electrons') > index)
            occ_beta = int(self.get('output/structure/dft/n_beta_electrons') > index)

            if self.get('output/structure/dft/beta_orbital_e') is None:
                orbital_energy = self.get('output/structure/dft/alpha_orbital_e')[index]
                print("#{}: \t Orbital energy = {:>10.5f} \t Occ. = {}".format(n,orbital_energy,occ_alpha+occ_beta))
            else:
                orbital_energy = [self.get('output/structure/dft/alpha_orbital_e')[index],self.get('output/structure/dft/beta_orbital_e')[index]]
                print("#{}: \t Orbital energies (alpha,beta) = {:>10.5f},{:>10.5f} \t Occ. = {},{}".format(n,orbital_energy[0],orbital_energy[1],occ_alpha,occ_beta))

    def visualize_MO(self, index, particle_size=0.5, show_bonds=True):
        """
            Visualize the MO identified by its index.

            **Arguments**

            index       index of the MO, as listed by print_MO()

            particle_size
                        size of the atoms for visualization, lower value if orbital is too small to see

            show_bonds  connect atoms or not

            **Notes**

            This function should always be accompanied with the following commands (in a separate cell)

            view[1].update_surface(isolevel=1, color='blue', opacity=.3)
            view[2].update_surface(isolevel=-1, color='red', opacity=.3)

            This makes sure that the bonding and non-bonding MO's are plotted and makes them transparent
        """
        n_MO = self.get('output/structure/dft/scf_density').shape[0]
        assert index >= 0 and index < n_MO
        assert len(self.get('output/structure/numbers')) < 50 # check whether structure does not become too large for interactive calculation of cube file

        # print orbital information
        occ_alpha = int(self.get('output/structure/dft/n_alpha_electrons') > index)
        occ_beta = int(self.get('output/structure/dft/n_beta_electrons') > index)

        if self.get('output/structure/dft/beta_orbital_e') is None:
            orbital_energy = self.get('output/structure/dft/alpha_orbital_e')[index]
            print("Orbital energy = {:>10.5f} \t Occ. = {}".format(orbital_energy,occ_alpha+occ_beta))
        else:
            orbital_energy = [self.get('output/structure/dft/alpha_orbital_e')[index],self.get('output/structure/dft/beta_orbital_e')[index]]
            print("Orbital energies (alpha,beta) = {:>10.5f},{:>10.5f} \t Occ. = {},{}".format(orbital_energy[0],orbital_energy[1],occ_alpha,occ_beta))

        # make cube file
        path = self.path+'_hdf5/'+self.name+'/input'
        out = subprocess.check_output(
                "ml load Gaussian/g16_E.01-intel-2019a;module use /apps/gent/CO7/haswell-ib/modules/all; cubegen 1 MO={} {}.fchk {}.cube".format(index+1,path,path),
                stderr=subprocess.STDOUT,
                universal_newlines=True,
                shell=True,
            )
        # visualize cube file
        try:
            import nglview
        except ImportError:
            raise ImportError("The animate_nma_mode() function requires the package nglview to be installed")

        atom_numbers = []
        atom_positions = []

        with open('{}.cube'.format(path),'r') as f:
            for i in range(2):
                f.readline()
            n_atoms = int(f.readline().split()[0][1:])
            for i in range(3):
                f.readline()
            for n in range(n_atoms):
                line = f.readline().split()
                atom_numbers.append(int(line[0]))
                atom_positions.append(np.array([float(m) for m in line[2:]])/angstrom)

        structure = Atoms(numbers=np.array(atom_numbers),positions=atom_positions)
        view = nglview.show_ase(structure)
        if not show_bonds:
            view.add_spacefill(radius_type='vdw', scale=0.5, radius=particle_size)
            view.remove_ball_and_stick()
        else:
            view.add_ball_and_stick()
        view.add_component('{}.cube'.format(path))
        view.add_component('{}.cube'.format(path))
        return view

    def read_NMA(self):
        """
            Reads the NMA output from the Gaussian .log file.

            Returns:
                    IR frequencies, intensities and corresponding eigenvectors (modes).
        """
        # Read number of atoms
        nrat = len(self.get('output/structure/numbers'))

        # Read IR frequencies and intensities from log file
        low_freqs = []
        freqs = []
        ints = []
        modes = [[] for i in range(nrat)]

        path = self.path+'_hdf5/'+self.name+'/input.log'
        with open(path,'r') as f:
            lines = f.readlines()

        # Assert normal termination
        assert "Normal termination of Gaussian" in lines[-1]

        # Find zero frequencies
        for n in range(len(lines)):
            line = lines[n]
            if 'Low frequencies' in line:
                low_freqs += [float(i) for i in line[20:].split()]
            if 'Frequencies --' in line:
                freqs += [float(i) for i in line[15:].split()]
            if 'IR Inten    --' in line:
                ints += [float(i) for i in line[15:].split()]
            if 'Atom  AN      X      Y      Z' in line:
                for m in range(nrat):
                    modes[m] += [float(i) for i in lines[n+m+1][10:].split()]

        nma_zeros = 3*nrat-len(freqs)
        freq_array = np.zeros(3*nrat)
        freq_array[:nma_zeros] = np.array(low_freqs[:nma_zeros])
        freq_array[nma_zeros:] = np.array(freqs)
        freqs = freq_array * (lightspeed/centimeter) # put into atomic units
        ints = np.array(ints)
        modes = np.array(modes).reshape(len(ints),nrat,3)

        return freqs,ints,modes

    def bsse_to_pandas(self):
        """
        Convert bsse output of all frames to a pandas Dataframe object.

        Returns:
            pandas.Dataframe: output as dataframe
        """
        assert 'counterpoise' in [k.lower() for k in self.input['settings'].keys()] # check if there was a bsse calculation
        tmp = {}
        with self.project_hdf5.open('output/structure/bsse') as hdf:
            for key in hdf.list_nodes():
                tmp[key] = hdf[key] if isinstance(hdf[key],np.ndarray) else [hdf[key]]
            df = pandas.DataFrame(tmp)
        return df


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


def write_input(input_dict,working_directory='.'):
    # Comments can be written with ! in Gaussian
    # Load dictionary
    lot          = input_dict['lot']
    basis_set    = input_dict['basis_set']
    spin_mult    = input_dict['spin_mult']  # 2S+1
    charge       = input_dict['charge']
    symbols      = input_dict['symbols']
    pos          = input_dict['pos']
    assert pos.shape[0] == len(symbols)

    # Optional elements
    if not input_dict['mem'] is None:
        mem = input_dict['mem'] + 'B' * (input_dict['mem'][-1]!='B') # check if string ends in bytes
        # convert pmem to mem
        cores = input_dict['cores']
        nmem = str(int(re.findall("\d+", mem)[0]) * cores)
        mem_unit = re.findall("[a-zA-Z]+", mem)[0]
        mem = nmem+mem_unit
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

    verbosity_dict = {'low': 't', 'normal': 'n', 'high': 'p'}
    if not input_dict['verbosity'] is None:
        verbosity  = input_dict['verbosity']
        if verbosity in verbosity_dict:
            verbosity = verbosity_dict[verbosity]
    else:
        verbosity='n'

    if 'Counterpoise' in settings.keys():
        if input_dict['bsse_idx'] is None or not len(input_dict['bsse_idx']) == len(pos) : # check if all elements are present for a BSSE calculation
            raise ValueError('The Counterpoise setting requires a valid bsse_idx array')
        # Check bsse idx (should start from 1 for Gaussian)
        input_dict['bsse_idx'] = [k - min(input_dict['bsse_idx']) + 1 for k in input_dict['bsse_idx']]
        # Check if it only contains conseqcutive numbers (sum of set should be n*(n+1)/2)
        assert sum(set(input_dict['bsse_idx'])) == (max(input_dict['bsse_idx'])*(max(input_dict['bsse_idx']) + 1))/2

    # Parse settings
    settings_string = ""
    for key,valuelst in settings.items():
        if not isinstance(valuelst, list):
            valuelst = [valuelst]
        option = key + "({}) ".format(",".join(valuelst))*(len(valuelst) > 0)
        settings_string += option

    # Write to file
    route_section = "#{} {}/{} {} {}\n\n".format(verbosity,lot,basis_set,jobtype,settings_string)
    with open(os.path.join(working_directory, 'input.com'), 'w') as f:
        f.write("%mem={}\n".format(mem))
        f.write("%chk=input.chk\n")
        f.write(route_section)
        f.write("{}\n\n".format(title))

        if not 'Counterpoise' in settings.keys():
            f.write("{} {}\n".format(charge,spin_mult))
            for n, p in enumerate(pos):
                f.write(" {}\t{: 1.6f}\t{: 1.6f}\t{: 1.6f}\n".format(symbols[n],p[0],p[1],p[2]))
            f.write('\n\n') # don't know whether this is still necessary in G16
        else:
            if isinstance(charge,list) and isinstance(spin_mult,list): # for BSSE it is possible to define charge and multiplicity for the fragments separately
                f.write(" ".join(["{},{}".format(charge[idx],spin_mult[idx]) for idx in range(int(settings['Counterpoise']))])) # first couple is for full system, then every fragment separately
            else:
                f.write("{} {}\n".format(charge,spin_mult))

            for n, p in enumerate(pos):
                f.write(" {}(Fragment={})\t{: 1.6f}\t{: 1.6f}\t{: 1.6f}\n".format(symbols[n],input_dict['bsse_idx'][n],p[0],p[1],p[2]))
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

    fchkdict['structure/numbers']     = fchk.fields.get('Atomic numbers')
    fchkdict['structure/masses']      = fchk.fields.get('Real atomic weights')*amu
    fchkdict['structure/charges']     = fchk.fields.get('Mulliken Charges')
    fchkdict['structure/dipole']      = fchk.fields.get('Dipole Moment')
    fchkdict['structure/dft/n_electrons']         = fchk.fields.get('Number of electrons')
    fchkdict['structure/dft/n_alpha_electrons']   = fchk.fields.get('Number of alpha electrons')
    fchkdict['structure/dft/n_beta_electrons']    = fchk.fields.get('Number of beta electrons')
    fchkdict['structure/dft/n_basis_functions']   = fchk.fields.get('Number of basis functions')

    # Orbital information
    fchkdict['structure/dft/alpha_orbital_e']     = fchk.fields.get('Alpha Orbital Energies')
    fchkdict['structure/dft/beta_orbital_e']      = fchk.fields.get('Beta Orbital Energies')

    # Densities
    fchkdict['structure/dft/scf_density']         = _triangle_to_dense(fchk.fields.get('Total SCF Density'))
    fchkdict['structure/dft/spin_scf_density']    = _triangle_to_dense(fchk.fields.get('Spin SCF Density'))

    if fchk.lot.upper() in ['MP2', 'MP3', 'CC', 'CI']:
        # only one of the lots should be present, hence using the same key
        fchkdict['structure/dft/post_scf_density']      = _triangle_to_dense(fchk.fields.get('Total {} Density'.format(fchk.lot)))
        fchkdict['structure/dft/post_spin_scf_density'] = _triangle_to_dense(fchk.fields.get('Spin {} Density'.format(fchk.lot)))

    # Specific job information
    if fchkdict['jobtype'] == 'fopt':
        if len(fchk.get_optimization_coordinates().shape) == 3:
            fchkdict['structure/positions']   = fchk.get_optimization_coordinates()[-1]/angstrom
        else:
            fchkdict['structure/positions']   = fchk.get_optimization_coordinates()/angstrom
        fchkdict['generic/positions']     = fchk.get_optimization_coordinates()/angstrom
        fchkdict['generic/energy_tot']    = fchk.get_optimization_energies()/electronvolt
        fchkdict['generic/forces']        = fchk.get_optimization_gradients()/(electronvolt/angstrom) * -1

    if fchkdict['jobtype'] == 'freq':
        fchkdict['structure/positions']   = fchk.fields.get('Current cartesian coordinates').reshape([1,-1, 3])/angstrom
        fchkdict['generic/positions']     = fchk.fields.get('Current cartesian coordinates').reshape([1,-1, 3])/angstrom
        fchkdict['generic/forces']        = fchk.fields.get('Cartesian Gradient').reshape([-1, 3])/(electronvolt/angstrom) *-1
        fchkdict['generic/hessian']       = fchk.get_hessian()/(electronvolt/angstrom**2)
        fchkdict['generic/energy_tot']    = fchk.fields.get('Total Energy')/electronvolt

    if fchkdict['jobtype'] == 'sp':
        fchkdict['structure/positions']   = fchk.fields.get('Current cartesian coordinates').reshape([1,-1, 3])/angstrom
        fchkdict['generic/positions']     = fchk.fields.get('Current cartesian coordinates').reshape([1,-1, 3])/angstrom
        fchkdict['generic/energy_tot']    = fchk.fields.get('Total Energy')/electronvolt

    return fchkdict


def get_bsse_array(line,it):
    numeric_const_pattern = '[-+]? (?: (?: \d* \. \d+ ) | (?: \d+ \.? ) )(?: [Ee] [+-]? \d+ ) ?'
    rx = re.compile(numeric_const_pattern, re.VERBOSE)

    cE_corr = float(rx.findall(line)[0]) * kcalmol/electronvolt
    line = next(it) # go to next line
    cE_raw = float(rx.findall(line)[0]) * kcalmol/electronvolt
    line = next(it) # go to next line
    sum_fragments = float(rx.findall(line)[0])/electronvolt
    line = next(it) # go to next line
    bsse_corr = float(rx.findall(line)[0])/electronvolt
    line = next(it) # go to next line
    E_tot_corr = float(rx.findall(line)[0])/electronvolt

    return E_tot_corr,bsse_corr,sum_fragments,cE_raw,cE_corr


def read_bsse(output_file,output_dict):
    # Check whether the route section contains the Counterpoise setting (if fchk module is update, route section can be loaded from dict)
    cp = False
    with open(output_file,'r') as f:
        line = f.readline()
        while line:
            if 'route' in line.lower():
                if 'counterpoise' in f.readline().lower(): # read next line
                    cp = True
                break
            line = f.readline()

    if cp:
        # the log file has the same path and name as the output file aside from the file extension
        log_file = output_file[:output_file.rfind('.')] + '.log'

        frames = 1 if isinstance(output_dict['generic/energy_tot'],float) else len(output_dict['generic/energy_tot'])

        output_dict['structure/bsse/energy_tot_corrected'] = np.zeros(frames)
        output_dict['structure/bsse/bsse_correction'] = np.zeros(frames)
        output_dict['structure/bsse/sum_of_fragments'] = np.zeros(frames)
        output_dict['structure/bsse/complexation_energy_raw'] = np.zeros(frames)
        output_dict['structure/bsse/complexation_energy_corrected'] = np.zeros(frames)

        it = _reverse_readline(log_file)
        line = next(it)
        for i in range(frames):
            found = False
            while not found:
                line = next(it)
                if 'complexation energy' in line:
                    E_tot_corr,bsse_corr,sum_fragments,cE_raw,cE_corr = get_bsse_array(line,it)
                    output_dict['structure/bsse/energy_tot_corrected'][i] = E_tot_corr
                    output_dict['structure/bsse/bsse_correction'][i] = bsse_corr
                    output_dict['structure/bsse/sum_of_fragments'][i] = sum_fragments
                    output_dict['structure/bsse/complexation_energy_raw'][i] = cE_raw
                    output_dict['structure/bsse/complexation_energy_corrected'][i] = cE_corr
                    found = True

        if frames==1:
            output_dict['structure/bsse/energy_tot_corrected'] = output_dict['structure/bsse/energy_tot_corrected'][0]
            output_dict['structure/bsse/bsse_correction'] = output_dict['structure/bsse/bsse_correction'][0]
            output_dict['structure/bsse/sum_of_fragments'] = output_dict['structure/bsse/sum_of_fragments'][0]
            output_dict['structure/bsse/complexation_energy_raw'] = output_dict['structure/bsse/complexation_energy_raw'][0]
            output_dict['structure/bsse/complexation_energy_corrected'] = output_dict['structure/bsse/complexation_energy_corrected'][0]
        else:
            # flip array sequence
            output_dict['structure/bsse/energy_tot_corrected'] = output_dict['structure/bsse/energy_tot_corrected'][::-1]
            output_dict['structure/bsse/bsse_correction'] = output_dict['structure/bsse/bsse_correction'][::-1]
            output_dict['structure/bsse/sum_of_fragments'] = output_dict['structure/bsse/sum_of_fragments'][::-1]
            output_dict['structure/bsse/complexation_energy_raw'] = output_dict['structure/bsse/complexation_energy_raw'][::-1]
            output_dict['structure/bsse/complexation_energy_corrected'] = output_dict['structure/bsse/complexation_energy_corrected'][::-1]


def read_EmpiricalDispersion(output_file,output_dict):
    # Get dispersion term from log file if it is there
    # dispersion term is not retrieved from gaussian output in fchk

    disp = None
    with open(output_file,'r') as f:
        while True:
            line = f.readline()
            if 'Route' in line:
                line = f.readline()
                if 'EmpiricalDispersion' in line:
                    idx = line.find('EmpiricalDispersion')
                    if 'GD3' in line[idx:]:
                        search_term = 'Grimme-D3 Dispersion energy='
                    else:
                        raise NotImplementedError
                else:
                    return
                break

    # the log file has the same path and name as the output file aside from the file extension
    log_file = output_file[:output_file.rfind('.')] + '.log'
    it = _reverse_readline(log_file)
    while True:
        line = next(it)
        if search_term in line:
            disp = float(line[38:-9])/electronvolt # could be changed when new search terms are implemented
            break

    output_dict['generic/energy_tot'] += disp


def collect_output(output_file):
    # Read output
    fchk = FCHKFile(output_file)

    # Translate to dict
    output_dict = fchk2dict(fchk)

    # Read BSSE output if it is present
    read_bsse(output_file,output_dict)

    # Correct energy if empirical dispersion contribution is present
    read_EmpiricalDispersion(output_file,output_dict)

    return output_dict


# function from theochem iodata
def _triangle_to_dense(triangle):
    """
    Convert a symmetric matrix in triangular storage to a dense square matrix.
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


def _reverse_readline(filename, buf_size=8192):
    """
    A generator that returns the lines of a file in reverse order

    https://stackoverflow.com/questions/2301789/read-a-file-in-reverse-order-using-python
    """
    with open(filename) as fh:
        segment = None
        offset = 0
        fh.seek(0, os.SEEK_END)
        file_size = remaining_size = fh.tell()
        while remaining_size > 0:
            offset = min(file_size, offset + buf_size)
            fh.seek(file_size - offset)
            buffer = fh.read(min(remaining_size, buf_size))
            remaining_size -= buf_size
            lines = buffer.split('\n')
            # The first line of the buffer is probably not a complete line so
            # we'll save it and append it to the last line of the next buffer
            # we read
            if segment is not None:
                # If the previous chunk starts right from the beginning of line
                # do not concat the segment to the last line of new chunk.
                # Instead, yield the segment first
                if buffer[-1] != '\n':
                    lines[-1] += segment
                else:
                    yield segment
            segment = lines[0]
            for index in range(len(lines) - 1, 0, -1):
                if lines[index]:
                    yield lines[index]
        # Don't yield None if the file was empty
        if segment is not None:
            yield segment
