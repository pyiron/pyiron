from pyiron_base import GenericParameters
from pyiron.atomistics.structure.atoms import Atoms
from pyiron.atomistics.job.atomistic import AtomisticGenericJob
from yaff import System, log as Yafflog, ForceField
Yafflog.set_level(Yafflog.silent)
from quickff import read_abinitio
from quickff.tools import set_ffatypes
from quickff.settings import key_checks
from molmod.units import *
from molmod.io.chk import load_chk, dump_chk
from molmod.periodic import periodic as pt

import os
import posixpath
import numpy as np


def write_chk(input_dict, working_directory='.'):
    # collect data and initialize Yaff system
    if 'cell' in input_dict.keys() and input_dict['cell'] is not None and np.all(np.array(input_dict['cell']) != np.zeros([3,3])):
        system = System(
          input_dict['numbers'], 
          input_dict['pos']*angstrom, 
          rvecs=np.array(input_dict['cell'])*angstrom, 
          ffatypes=input_dict['ffatypes_man'], 
          ffatype_ids=input_dict['ffatype_ids_man']
        )
    else:
        system = System(
            input_dict['numbers'],
            input_dict['pos']*angstrom,
            ffatypes=input_dict['ffatypes_man'],
            ffatype_ids=input_dict['ffatype_ids_man']
        )
    # determine masses, bonds and ffaypes from ffatype_rules
    system.detect_bonds()
    system.set_standard_masses()
    # write dictionnairy to MolMod CHK file
    system.to_file(posixpath.join(working_directory, 'input.chk'))
    # Reload input.chk as dictionairy and add AI input data
    d = load_chk(posixpath.join(working_directory, 'input.chk'))

    assert isinstance(input_dict['aiener'], float), "AI energy not defined in input, use job.read_abintio(...)"
    assert isinstance(input_dict['aigrad'], np.ndarray), "AI gradient not defined in input, use job.read_abintio(...)"
    assert isinstance(input_dict['aihess'], np.ndarray), "AI hessian not defined in input, use job.read_abintio(...)"
    d['energy'] = input_dict['aiener']
    d['grad'] = input_dict['aigrad']
    d['hess'] = input_dict['aihess']
    dump_chk(posixpath.join(working_directory,'input.chk'), d)


def write_pars(pars, fn, working_directory='.'):
    with open(posixpath.join(working_directory, fn), 'w') as f:
        for line in pars:
            f.write(line)


def write_config(input_dict, working_directory='.'):
    with open(posixpath.join(working_directory, 'config.txt'), 'w') as f:
        for key in key_checks.keys():
            if key in input_dict.keys():
                value_int = str(input_dict[key])
                if key == 'ffatypes': assert value_int == 'None'
                print('%s:   %s' %(key+' '*(30-len(key)), value_int), file=f)


def collect_output(fn_pars, fn_sys):
    # this routine reads the output parameter file containing the covalent pars
    output_dict = {
        'generic/bond': [],
        'generic/bend': [],
        'generic/torsion': [],
        'generic/oopdist': [],
        'generic/cross': []
    }
    kinds = ['bond', 'bend', 'torsion', 'oopdist', 'cross']
    with open(fn_pars, 'r') as f:
        for line in f.readlines():
            for key in kinds:
                if key in line.lower():
                    output_dict['generic/%s' %key].append(line)
    system = System.from_file(fn_sys)
    output_dict['system/numbers'] = system.numbers
    output_dict['system/pos'] = system.pos/angstrom
    if system.cell is not None:
        output_dict['system/rvecs'] = system.cell.rvecs/angstrom
    output_dict['system/bonds'] = system.bonds
    output_dict['system/ffatypes'] = np.asarray(system.ffatypes,'S22')
    output_dict['system/ffatype_ids'] = system.ffatype_ids
    return output_dict


class QuickFFInput(GenericParameters):
    def __init__(self, input_file_name=None):
        super(QuickFFInput, self).__init__(input_file_name=input_file_name,table_name="input_inp",comment_char="#")

    def load_default(self):
        """
        Loading the default settings for the input file.
        """
        input_str = """\
fn_yaff pars_cov.txt
fn_charmm22_prm None
fn_charmm22_psf None
fn_sys system.chk
plot_traj None
xyz_traj False
fn_traj None
log_level high
log_file quickff.log
program_mode DeriveFF
only_traj PT_ALL
ffatypes None #Define atom types using the built-in routine in QuickFF (see documentation)
ei None
ei_rcut None #default is 20 (periodic) or 50 (non-per) A
vdw None
vdw_rcut 37.79452267842504
covres None
excl_bonds None
excl_bends None
excl_dihs None
excl_oopds None
do_hess_mass_weighting True
do_hess_negfreq_proj False
do_cross_svd True
pert_traj_tol 1e-3
pert_traj_energy_noise None
cross_svd_rcond 1e-8
do_bonds True
do_bends True
do_dihedrals True
do_oops True
do_cross_ASS True
do_cross_ASA True
do_cross_DSS False
do_cross_DSD False
do_cross_DAA False
do_cross_DAD False
consistent_cross_rvs True
remove_dysfunctional_cross True
bond_term BondHarm
bend_term BendAHarm
do_squarebend True
do_bendclin True
do_sqoopdist_to_oopdist True
"""
        self.load_string(input_str)


class QuickFF(AtomisticGenericJob):
    def __init__(self, project, job_name):
        super(QuickFF, self).__init__(project, job_name)
        self.__name__ = "QuickFF"
        self._executable_activate(enforce=True)
        self.input = QuickFFInput()
        self.ffatypes = None
        self.ffatype_ids = None
        self.aiener = None
        self.aigrad = None
        self.aihess = None
        self.fn_ei = None
        self.fn_vdw = None

    def read_abinitio(self, fn):
        numbers, coords, energy, grad, hess, masses, rvecs, pbc = read_abinitio(fn)
        coords /= angstrom
        if rvecs is not None:
            rvecs /= angstrom
        self.structure = Atoms(numbers=numbers, positions=coords, cell=rvecs, pbc=True)
        self.aiener = energy
        self.aigrad = grad
        self.aihess = hess

    def detect_ffatypes(self, ffatypes=None, ffatype_rules=None, ffatype_level=None):
        """
        Define atom types by explicitely giving them through the
        ffatypes keyword, defining atype rules using the ATSELECT
        language implemented in Yaff (see the Yaff documentation at
        http://molmod.github.io/yaff/ug_atselect.html) or by specifying
        the ffatype_level employing the built-in routine in QuickFF.
        """
        numbers = np.array([pt[symbol].number for symbol in self.structure.get_chemical_symbols()])
        if self.structure.cell is not None and np.all(np.array(self.structure.cell) != np.zeros([3,3])):
            system = System(numbers, self.structure.positions.copy()*angstrom, rvecs=np.array(self.structure.cell)*angstrom)
        else:
            system = System(numbers, self.structure.positions.copy()*angstrom)
        system.detect_bonds()

        if not sum([ffatypes is None, ffatype_rules is None, ffatype_level is None]) == 2:
            raise IOError('Exactly one of ffatypes, ffatype_rules and ffatype_level should be defined')

        if ffatypes is not None:
            assert ffatype_rules is None, 'ffatypes and ffatype_rules cannot be defined both'
            system.ffatypes = ffatypes
            system.ffatype_ids = None
            system._init_derived_ffatypes()
        if ffatype_rules is not None:
            system.detect_ffatypes(ffatype_rules)
        if ffatype_level is not None:
            set_ffatypes(system, ffatype_level)

        self.ffatypes = system.ffatypes.copy()
        self.ffatype_ids = system.ffatype_ids.copy()

    def set_ei(self, fn):
        self.input['ei'] = fn.split('/')[-1]
        self.fn_ei = fn

    def set_vdw(self, fn):
        self.input['vdw'] = fn.split('/')[-1]
        self.fn_vdw = fn

    def write_input(self):
        # load system related input
        input_dict = {
            'symbols': self.structure.get_chemical_symbols(),
            'numbers': np.array([pt[symbol].number for symbol in self.structure.get_chemical_symbols()]),
            'ffatypes_man': self.ffatypes,
            'ffatype_ids_man': self.ffatype_ids,
            'pos': self.structure.positions,
            'aiener': self.aiener,
            'aigrad': self.aigrad,
            'aihess': self.aihess,
        }
        for k in self.input._dataset["Parameter"]:
            input_dict[k] = self.input[k]
        input_dict['cell'] = None
        if self.structure.cell is not None:
            input_dict['cell'] = self.structure.get_cell()
        # load all input settings from self.input
        for k, v in self.input._dataset.items():
            input_dict[k] = v
        # write input chk file
        write_chk(input_dict, working_directory=self.working_directory)
        # write nonbonded pars and config input files
        if self.fn_ei is not None:
            assert self.input['ei'] is not None
            os.system('cp %s %s/%s' %(self.fn_ei , self.working_directory, self.input['ei']))
        if self.fn_vdw is not None:
            assert self.input['vdw'] is not None
            os.system('cp %s %s/%s' %(self.fn_vdw, self.working_directory, self.input['vdw']))
        write_config(input_dict,working_directory=self.working_directory)

    def collect_output(self):
        output_dict = collect_output(
            posixpath.join(self.working_directory, self.input['fn_yaff']),
            posixpath.join(self.working_directory, self.input['fn_sys'])
        )
        with self.project_hdf5.open("output") as hdf5_output:
            for k, v in output_dict.items():
                hdf5_output[k] = v

    def to_hdf(self, hdf=None, group_name=None):
        super(QuickFF, self).to_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf5_input:
            self.structure.to_hdf(hdf5_input)
            self.input.to_hdf(hdf5_input)

    def from_hdf(self, hdf=None, group_name=None):
        super(QuickFF, self).from_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf5_input:
            self.input.from_hdf(hdf5_input)
            self.structure = Atoms().from_hdf(hdf5_input)

    def get_structure(self, iteration_step=-1, wrap_atoms=True):
        """
        Overwrite the get_structure routine from AtomisticGenericJob because we want to avoid
        defining a unit cell when one does not exist
        """
        raise NotImplementedError

    def log(self):
        with open(posixpath.join(self.working_directory, 'quickff.log')) as f:
            print(f.read())

    def get_yaff_system(self):
        system = System.from_file(posixpath.join(self.working_directory, self.input['fn_sys']))
        return system

    def get_yaff_ff(self, rcut=15*angstrom, alpha_scale=3.2, gcut_scale=1.5, smooth_ei=True):
        system = self.get_yaff_system()
        fn_pars = posixpath.join(self.working_directory, self.input['fn_yaff'])
        if not os.path.isfile(fn_pars):
            raise IOError('No pars.txt file find in job working directory. Have you already run the job?')
        ff = ForceField.generate(
            system, fn_pars, rcut=rcut, alpha_scale=alpha_scale,
            gcut_scale=gcut_scale, smooth_ei=smooth_ei
        )
        return ff
