from __future__ import print_function
# lammps.py (2011/03/29)
# An ASE calculator for the LAMMPS classical MD code available from
#       http://lammps.sandia.gov/
# The environment variable LAMMPS_COMMAND must be defined to point to the LAMMPS binary.
#
# Copyright (C) 2009 - 2011 Joerg Meyer, joerg.meyer@ch.tum.de
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this file; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
# or see <http://www.gnu.org/licenses/>.


import os
import shutil
import shlex
from subprocess import Popen, PIPE
from threading import Thread
from re import compile as re_compile, IGNORECASE
from tempfile import mkdtemp, NamedTemporaryFile, mktemp as uns_mktemp
import numpy as np
import decimal as dec
from ase import Atoms
from ase.parallel import paropen
from ase.units import GPa

__all__ = ['LAMMPS', 'write_lammps_data']

# "End mark" used to indicate that the calculation is done
CALCULATION_END_MARK = '__end_of_ase_invoked_calculation__'


class LAMMPS:
    def __init__(self, label='lammps', tmp_dir=None, parameters={},
                 specorder=None, files=[], always_triclinic=False,
                 keep_alive=True, keep_tmp_files=False,
                 no_data_file=False):
        """The LAMMPS calculators object

        files: list
            Short explanation XXX
        parameters: dict
            Short explanation XXX
        specorder: list
            Short explanation XXX
        keep_tmp_files: bool
            Retain any temporary files created. Mostly useful for debugging.
        tmp_dir: str
            path/dirname (default None -> create automatically).
            Explicitly control where the calculator object should create
            its files. Using this option implies 'keep_tmp_files'
        no_data_file: bool
            Controls whether an explicit data file will be used for feeding
            atom coordinates into lammps. Enable it to lessen the pressure on
            the (tmp) file system. THIS OPTION MIGHT BE UNRELIABLE FOR CERTAIN
            CORNER CASES (however, if it fails, you will notice...).
        keep_alive: bool
            When using LAMMPS as a spawned subprocess, keep the subprocess
            alive (but idling when unused) along with the calculator object.
        always_triclinic: bool
            Force use of a triclinic cell in LAMMPS, even if the cell is
            a perfect parallelepiped.
        """

        self.label = label
        self.parameters = parameters
        self.specorder = specorder
        self.files = files
        self.always_triclinic = always_triclinic
        self.calls = 0
        self.forces = None
        self.keep_alive = keep_alive
        self.keep_tmp_files = keep_tmp_files
        self.no_data_file = no_data_file
        if tmp_dir is not None:
            # If tmp_dir is pointing somewhere, don't remove stuff!
            self.keep_tmp_files = True
        self._lmp_handle = None  # To handle the lmp process

        # read_log depends on that the first (three) thermo_style custom args
        # can be capitilized and matched against the log output. I.e.
        # don't use e.g. 'ke' or 'cpu' which are labeled KinEng and CPU.
        self._custom_thermo_args = ['step', 'temp', 'press', 'cpu',
                                    'pxx', 'pyy', 'pzz', 'pxy', 'pxz', 'pyz',
                                    'ke', 'pe', 'etotal',
                                    'vol', 'lx', 'ly', 'lz', 'atoms']
        self._custom_thermo_mark = ' '.join([x.capitalize() for x in
                                             self._custom_thermo_args[0:3]])

        # Match something which can be converted to a float
        f_re = r'([+-]?(?:(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?|nan|inf))'
        n = len(self._custom_thermo_args)
        # Create a re matching exactly N white space separated floatish things
        self._custom_thermo_re = re_compile(r'^\s*' + r'\s+'.join([f_re] * n) + r'\s*$',
                                            flags=IGNORECASE)
        # thermo_content contains data "written by" thermo_style.
        # It is a list of dictionaries, each dict (one for each line
        # printed by thermo_style) contains a mapping between each
        # custom_thermo_args-argument and the corresponding
        # value as printed by lammps. thermo_content will be
        # re-populated by the read_log method.
        self.thermo_content = []

        if tmp_dir is None:
            self.tmp_dir = mkdtemp(prefix='LAMMPS-')
        else:
            self.tmp_dir = os.path.realpath(tmp_dir)
            if not os.path.isdir(self.tmp_dir):
                os.mkdir(self.tmp_dir, 0o755)

        for f in files:
            shutil.copy(f, os.path.join(self.tmp_dir, os.path.basename(f)))

    def clean(self, force=False):

        self._lmp_end()

        if not self.keep_tmp_files:
            shutil.rmtree(self.tmp_dir)

    def get_potential_energy(self, atoms):
        self.update(atoms)
        return self.thermo_content[-1]['pe']

    def get_forces(self, atoms):
        self.update(atoms)
        return self.forces.copy()

    def get_stress(self, atoms):
        self.update(atoms)
        tc = self.thermo_content[-1]
        # 1 bar (used by lammps for metal units) = 1e-4 GPa
        return np.array([tc[i] for i in ('pxx', 'pyy', 'pzz',
                                         'pyz', 'pxz', 'pxy')]) * (-1e-4 * GPa)

    def update(self, atoms):
        if not hasattr(self, 'atoms') or self.atoms != atoms:
            self.calculate(atoms)

    def calculate(self, atoms):
        self.atoms = atoms.copy()
        pbc = self.atoms.get_pbc()
        if all(pbc):
            cell = self.atoms.get_cell()
        elif not any(pbc):
            # large enough cell for non-periodic calculation -
            # LAMMPS shrink-wraps automatically via input command
            #       "periodic s s s"
            # below
            cell = 2 * np.max(np.abs(self.atoms.get_positions())) * np.eye(3)
        else:
            print("WARNING: semi-periodic ASE cell detected -")
            print("         translation to proper LAMMPS input cell might fail")
            cell = self.atoms.get_cell()
        self.prism = prism(cell)
        self.run()

    def _lmp_alive(self):
        # Return True if this calculator is currently handling a running lammps process
        return self._lmp_handle and not isinstance(self._lmp_handle.poll(), int)

    def _lmp_end(self):
        # Close lammps input and wait for lammps to end. Return process return value
        if self._lmp_alive():
            self._lmp_handle.stdin.close()
            return self._lmp_handle.wait()

    def run(self):
        """Method which explicitely runs LAMMPS."""

        self.calls += 1

        # set LAMMPS command from environment variable
        if 'LAMMPS_COMMAND' in os.environ:
            lammps_cmd_line = shlex.split(os.environ['LAMMPS_COMMAND'])
            if len(lammps_cmd_line) == 0:
                self.clean()
                raise RuntimeError('The LAMMPS_COMMAND environment variable '
                                   'must not be empty')
            # want always an absolute path to LAMMPS binary when calling from self.dir
            lammps_cmd_line[0] = os.path.abspath(lammps_cmd_line[0])

        else:
            self.clean()
            raise RuntimeError('Please set LAMMPS_COMMAND environment variable')
        if 'LAMMPS_OPTIONS' in os.environ:
            lammps_options = shlex.split(os.environ['LAMMPS_OPTIONS'])
        else:
            lammps_options = shlex.split('-echo log -screen none')

        # change into subdirectory for LAMMPS calculations
        cwd = os.getcwd()
        os.chdir(self.tmp_dir)

        # setup file names for LAMMPS calculation
        label = '{0}{1:>06}'.format(self.label, self.calls)
        lammps_in = uns_mktemp(prefix='in_' + label, dir=self.tmp_dir)
        lammps_log = uns_mktemp(prefix='log_' + label, dir=self.tmp_dir)
        lammps_trj_fd = NamedTemporaryFile(prefix='trj_' + label, dir=self.tmp_dir,
                                           delete=(not self.keep_tmp_files))
        lammps_trj = lammps_trj_fd.name
        if self.no_data_file:
            lammps_data = None
        else:
            lammps_data_fd = NamedTemporaryFile(prefix='data_' + label, dir=self.tmp_dir,
                                                delete=(not self.keep_tmp_files))
            self.write_lammps_data(lammps_data=lammps_data_fd)
            lammps_data = lammps_data_fd.name
            lammps_data_fd.flush()

        # see to it that LAMMPS is started
        if not self._lmp_alive():
            # Attempt to (re)start lammps
            self._lmp_handle = Popen(lammps_cmd_line + lammps_options + ['-log', '/dev/stdout'],
                                     stdin=PIPE, stdout=PIPE)
        lmp_handle = self._lmp_handle

        # Create thread reading lammps stdout (for reference, if requested,
        # also create lammps_log, although it is never used)
        if self.keep_tmp_files:
            lammps_log_fd = open(lammps_log, 'wb')
            fd = special_tee(lmp_handle.stdout, lammps_log_fd)
        else:
            fd = lmp_handle.stdout
        thr_read_log = Thread(target=self.read_lammps_log, args=(fd,))
        thr_read_log.start()

        # write LAMMPS input (for reference, also create the file lammps_in,
        # although it is never used)
        if self.keep_tmp_files:
            lammps_in_fd = open(lammps_in, 'wb')
            fd = special_tee(lmp_handle.stdin, lammps_in_fd)
        else:
            fd = lmp_handle.stdin
        self.write_lammps_in(lammps_in=fd, lammps_trj=lammps_trj, lammps_data=lammps_data)

        if self.keep_tmp_files:
            lammps_in_fd.close()

        # Wait for log output to be read (i.e., for LAMMPS to finish)
        # and close the log file if there is one
        thr_read_log.join()
        if self.keep_tmp_files:
            lammps_log_fd.close()

        if not self.keep_alive:
            self._lmp_end()

        exitcode = lmp_handle.poll()
        if exitcode and exitcode != 0:
            cwd = os.getcwd()
            raise RuntimeError('LAMMPS exited in {0} with exit code: {0}.'
                               ''.format(cwd, exitcode))

        # A few sanity checks
        if len(self.thermo_content) == 0:
            raise RuntimeError('Failed to retrieve any thermo_style-output')
        if int(self.thermo_content[-1]['atoms']) != len(self.atoms):
            # This obviously shouldn't happen, but if prism.fold_...() fails, it could
            raise RuntimeError('Atoms have gone missing')

        self.read_lammps_trj(lammps_trj=lammps_trj)
        lammps_trj_fd.close()
        if not self.no_data_file:
            lammps_data_fd.close()

        os.chdir(cwd)

    def write_lammps_data(self, lammps_data=None):
        """Method which writes a LAMMPS data file with atomic structure."""
        if (lammps_data == None):
            lammps_data = 'data.' + self.label
        write_lammps_data(lammps_data, self.atoms, self.specorder,
                          force_skew=self.always_triclinic, prismobj=self.prism)

    def write_lammps_in(self, lammps_in=None, lammps_trj=None, lammps_data=None):
        """Method which writes a LAMMPS in file with run parameters and settings."""

        if isinstance(lammps_in, str):
            f = paropen(lammps_in, 'wb')
            close_in_file = True
        else:
            # Expect lammps_in to be a file-like object
            f = lammps_in
            close_in_file = False

        if self.keep_tmp_files:
            f.write('# (written by ASE)\n'.encode('utf-8'))

        # Write variables
        f.write(('clear\n'
                 'variable dump_file string "{0}"\n'
                 'variable data_file string "{1}"\n'
                 ).format(lammps_trj, lammps_data).encode('utf-8'))

        parameters = self.parameters
        pbc = self.atoms.get_pbc()
        f.write('units metal \n'.encode('utf-8'))
        if ('boundary' in parameters):
            f.write('boundary {0} \n'.format(parameters['boundary']
                                             ).encode('utf-8'))
        else:
            f.write('boundary {0} {1} {2} \n'.format(
                *tuple('sp'[x] for x in pbc)).encode('utf-8'))
        f.write('atom_modify sort 0 0.0 \n'.encode('utf-8'))
        for key in ('neighbor', 'newton'):
            if key in parameters:
                f.write('{0} {1} \n'.format(key, parameters[key]
                                            ).encode('utf-8'))
        f.write('\n'.encode('utf-8'))

        # If self.no_lammps_data,
        # write the simulation box and the atoms
        if self.no_data_file:
            if self.keep_tmp_files:
                f.write('## Original ase cell\n'.encode('utf-8'))
                f.write(''.join(['# {0:.16} {1:.16} {2:.16}\n'.format(*x)
                                 for x in self.atoms.get_cell()]
                                ).encode('utf-8'))

            p = self.prism
            f.write('lattice sc 1.0\n'.encode('utf-8'))
            xhi, yhi, zhi, xy, xz, yz = p.get_lammps_prism_str()
            if self.always_triclinic or p.is_skewed():
                f.write('region asecell prism 0.0 {0} 0.0 {1} 0.0 {2} '.format(
                    xhi, yhi, zhi).encode('utf-8'))
                f.write('{0} {1} {2} side in units box\n'.format(
                    xy, xz, yz).encode('utf-8'))
            else:
                f.write(('region asecell block 0.0 {0} 0.0 {1} 0.0 {2} '
                         'side in units box\n').format(xhi, yhi, zhi
                                                       ).encode('utf-8'))

            symbols = self.atoms.get_chemical_symbols()
            if self.specorder is None:
                # By default, atom types in alphabetic order
                species = sorted(set(symbols))
            else:
                # By request, specific atom type ordering
                species = self.specorder

            n_atom_types = len(species)
            species_i = dict([(s, i + 1) for i, s in enumerate(species)])

            f.write('create_box {0} asecell\n'.format(n_atom_types
                                                      ).encode('utf-8'))
            for s, pos in zip(symbols, self.atoms.get_positions()):
                if self.keep_tmp_files:
                    f.write('# atom pos in ase cell: {0:.16} {1:.16} {2:.16}\n'
                            ''.format(*tuple(pos)).encode('utf-8'))
                f.write('create_atoms {0} single {1} {2} {3} units box\n'.format(
                    *((species_i[s],) + p.pos_to_lammps_fold_str(pos))
                ).encode('utf-8'))


        # if NOT self.no_lammps_data, then simply refer to the data-file
        else:
            f.write('read_data "{0}"\n'.format(lammps_data).encode('utf-8'))
            #f.write('read_data $data_file\n'.encode('utf-8'))

        # Write interaction stuff
        f.write('\n### interactions \n'.encode('utf-8'))
        if (('pair_style' in parameters) and ('pair_coeff' in parameters)):
            pair_style = parameters['pair_style']
            f.write('pair_style {0} \n'.format(pair_style).encode('utf-8'))
            for pair_coeff in parameters['pair_coeff']:
                f.write('pair_coeff {0} \n'.format(pair_coeff).encode('utf-8'))
            if 'mass' in parameters:
                for mass in parameters['mass']:
                    f.write('mass {0} \n'.format(mass).encode('utf-8'))
        else:
            # simple default parameters
            # that should always make the LAMMPS calculation run
            f.write('pair_style lj/cut 2.5 \n'
                    'pair_coeff * * 1 1 \n'
                    'mass * 1.0 \n'.encode('utf-8'))

        f.write('\n### run\n'
                'fix fix_nve all nve\n'
                'dump dump_all all custom 1 "{0}" id type x y z vx vy vz fx fy fz\n'
                ''.format(lammps_trj).encode('utf-8'))

        f.write('thermo_style custom {0}\n'
                'thermo_modify flush yes\n'
                'thermo 1\n'.format(' '.join(self._custom_thermo_args)
                                    ).encode('utf-8'))

        if 'minimize' in parameters:
            f.write('minimize {0}\n'.format(parameters['minimize']).encode('utf-8'))
        if 'run' in parameters:
            f.write('run {0}\n'.format(parameters['run']).encode('utf-8'))
        if not (('minimize' in parameters) or ('run' in parameters)):
            f.write('run 0\n'.encode('utf-8'))

        f.write('print "{0}"\n'.format(CALCULATION_END_MARK).encode('utf-8'))
        f.write('log /dev/stdout\n'.encode('utf-8'))  # Force LAMMPS to flush log

        f.flush()
        if close_in_file:
            f.close()

    def read_lammps_log(self, lammps_log=None, PotEng_first=False):
        """Method which reads a LAMMPS output log file."""

        if (lammps_log == None):
            lammps_log = self.label + '.log'

        if isinstance(lammps_log, str):
            f = paropen(lammps_log, 'wb')
            close_log_file = True
        else:
            # Expect lammps_in to be a file-like object
            f = lammps_log
            close_log_file = False

        thermo_content = []
        line = f.readline().decode('utf-8')
        while line and line.strip() != CALCULATION_END_MARK:
            # get thermo output
            if line.startswith(self._custom_thermo_mark):
                m = True
                while m:
                    line = f.readline().decode('utf-8')
                    m = self._custom_thermo_re.match(line)
                    if m:
                        # create a dictionary between each of the thermo_style args
                        # and it's corresponding value
                        thermo_content.append(dict(zip(self._custom_thermo_args,
                                                       map(float, m.groups()))))
            else:
                line = f.readline().decode('utf-8')

        if close_log_file:
            f.close()

        self.thermo_content = thermo_content

    def read_lammps_trj(self, lammps_trj=None, set_atoms=False):
        """Method which reads a LAMMPS dump file."""
        if (lammps_trj == None):
            lammps_trj = self.label + '.lammpstrj'

        f = paropen(lammps_trj, 'r')
        while True:
            line = f.readline()

            if not line:
                break

            # TODO: extend to proper dealing with multiple steps in one trajectory file
            if 'ITEM: TIMESTEP' in line:
                n_atoms = 0
                lo = [];
                hi = [];
                tilt = []
                id = [];
                type = []
                positions = [];
                velocities = [];
                forces = []

            if 'ITEM: NUMBER OF ATOMS' in line:
                line = f.readline()
                n_atoms = int(line.split()[0])

            if 'ITEM: BOX BOUNDS' in line:
                # save labels behind "ITEM: BOX BOUNDS" in triclinic case (>=lammps-7Jul09)
                tilt_items = line.split()[3:]
                for i in range(3):
                    line = f.readline()
                    fields = line.split()
                    lo.append(float(fields[0]))
                    hi.append(float(fields[1]))
                    if (len(fields) >= 3):
                        tilt.append(float(fields[2]))

            if 'ITEM: ATOMS' in line:
                # (reliably) identify values by labels behind "ITEM: ATOMS" - requires >=lammps-7Jul09
                # create corresponding index dictionary before iterating over atoms to (hopefully) speed up lookups...
                atom_attributes = {}
                for (i, x) in enumerate(line.split()[2:]):
                    atom_attributes[x] = i
                for n in range(n_atoms):
                    line = f.readline()
                    fields = line.split()
                    id.append(int(fields[atom_attributes['id']]))
                    type.append(int(fields[atom_attributes['type']]))
                    positions.append([float(fields[atom_attributes[x]]) for x in ['x', 'y', 'z']])
                    velocities.append([float(fields[atom_attributes[x]]) for x in ['vx', 'vy', 'vz']])
                    forces.append([float(fields[atom_attributes[x]]) for x in ['fx', 'fy', 'fz']])
        f.close()

        # determine cell tilt (triclinic case!)
        if (len(tilt) >= 3):
            # for >=lammps-7Jul09 use labels behind "ITEM: BOX BOUNDS" to assign tilt (vector) elements ...
            if (len(tilt_items) >= 3):
                xy = tilt[tilt_items.index('xy')]
                xz = tilt[tilt_items.index('xz')]
                yz = tilt[tilt_items.index('yz')]
            # ... otherwise assume default order in 3rd column (if the latter was present)
            else:
                xy = tilt[0]
                xz = tilt[1]
                yz = tilt[2]
        else:
            xy = xz = yz = 0
        xhilo = (hi[0] - lo[0]) - xy - xz
        yhilo = (hi[1] - lo[1]) - yz
        zhilo = (hi[2] - lo[2])

        # The simulation box bounds are included in each snapshot and if the box is triclinic (non-orthogonal),
        # then the tilt factors are also printed; see the region prism command for a description of tilt factors.
        # For triclinic boxes the box bounds themselves (first 2 quantities on each line) are a true "bounding box"
        # around the simulation domain, which means they include the effect of any tilt.
        # [ http://lammps.sandia.gov/doc/dump.html , lammps-7Jul09 ]
        #
        # This *should* extract the lattice vectors that LAMMPS uses from the true "bounding box" printed in the dump file
        # It might fail in some cases (negative tilts?!) due to the MIN / MAX construction of these box corners:
        #
        #       void Domain::set_global_box()
        #       [...]
        #         if (triclinic) {
        #           [...]
        #           boxlo_bound[0] = MIN(boxlo[0],boxlo[0]+xy);
        #           boxlo_bound[0] = MIN(boxlo_bound[0],boxlo_bound[0]+xz);
        #           boxlo_bound[1] = MIN(boxlo[1],boxlo[1]+yz);
        #           boxlo_bound[2] = boxlo[2];
        #
        #           boxhi_bound[0] = MAX(boxhi[0],boxhi[0]+xy);
        #           boxhi_bound[0] = MAX(boxhi_bound[0],boxhi_bound[0]+xz);
        #           boxhi_bound[1] = MAX(boxhi[1],boxhi[1]+yz);
        #           boxhi_bound[2] = boxhi[2];
        #         }
        # [ lammps-7Jul09/src/domain.cpp ]
        #
        cell = [[xhilo, 0, 0], [xy, yhilo, 0], [xz, yz, zhilo]]

        # assume that LAMMPS does not reorder atoms internally
        cell_atoms = np.array(cell)
        type_atoms = np.array(type)

        if self.atoms:
            cell_atoms = self.atoms.get_cell()

            # BEWARE: reconstructing the rotation from the LAMMPS output trajectory file
            #         fails in case of shrink wrapping for a non-periodic direction
            # -> hence rather obtain rotation from prism object used to generate the LAMMPS input
            # rotation_lammps2ase = np.dot(np.linalg.inv(np.array(cell)), cell_atoms)
            rotation_lammps2ase = np.linalg.inv(self.prism.R)

            type_atoms = self.atoms.get_atomic_numbers()
            positions_atoms = np.array([np.dot(np.array(r), rotation_lammps2ase) for r in positions])
            # velocities_atoms = np.array( [np.dot(np.array(v), rotation_lammps2ase) for v in velocities] )
            forces_atoms = np.array([np.dot(np.array(f), rotation_lammps2ase) for f in forces])

        if (set_atoms):
            # assume periodic boundary conditions here (like also below in write_lammps)
            self.atoms = Atoms(type_atoms, positions=positions_atoms, cell=cell_atoms)

        self.forces = forces_atoms


class special_tee:
    """A special purpose, with limited applicability, tee-like thing.

    A subset of stuff read from, or written to, orig_fd,
    is also written to out_fd.
    It is used by the lammps calculator for creating file-logs of stuff read from,
    or written to, stdin and stdout, respectively.
    """

    def __init__(self, orig_fd, out_fd):
        self._orig_fd = orig_fd
        self._out_fd = out_fd
        self.name = orig_fd.name

    def write(self, data):
        self._orig_fd.write(data)
        self._out_fd.write(data)
        self.flush()

    def read(self, *args, **kwargs):
        data = self._orig_fd.read(*args, **kwargs)
        self._out_fd.write(data)
        return data

    def readline(self, *args, **kwargs):
        data = self._orig_fd.readline(*args, **kwargs)
        self._out_fd.write(data)
        return data

    def readlines(self, *args, **kwargs):
        data = self._orig_fd.readlines(*args, **kwargs)
        self._out_fd.write(''.join(data))
        return data

    def flush(self):
        self._orig_fd.flush()
        self._out_fd.flush()


class prism:
    def __init__(self, cell, pbc=(True, True, True), digits=10):
        """Create a lammps-style triclinic prism object from a cell

        The main purpose of the prism-object is to create suitable
        string representations of prism limits and atom positions
        within the prism.
        When creating the object, the digits parameter (default set to 10)
        specify the precission to use.
        lammps is picky about stuff being within semi-open intervals,
        e.g. for atom positions (when using create_atom in the in-file),
        x must be within [xlo, xhi).
        """
        a, b, c = cell
        an, bn, cn = [np.linalg.norm(v) for v in cell]

        alpha = np.arccos(np.dot(b, c) / (bn * cn))
        beta = np.arccos(np.dot(a, c) / (an * cn))
        gamma = np.arccos(np.dot(a, b) / (an * bn))

        xhi = an
        xyp = np.cos(gamma) * bn
        yhi = np.sin(gamma) * bn
        xzp = np.cos(beta) * cn
        yzp = (bn * cn * np.cos(alpha) - xyp * xzp) / yhi
        zhi = np.sqrt(cn ** 2 - xzp ** 2 - yzp ** 2)

        # Set precision
        self.car_prec = dec.Decimal('10.0') ** \
                        int(np.floor(np.log10(max((xhi, yhi, zhi)))) - digits)
        self.dir_prec = dec.Decimal('10.0') ** (-digits)
        self.acc = float(self.car_prec)
        self.eps = np.finfo(xhi).eps

        # For rotating positions from ase to lammps
        Apre = np.array(((xhi, 0, 0),
                         (xyp, yhi, 0),
                         (xzp, yzp, zhi)))
        self.R = np.dot(np.linalg.inv(cell), Apre)

        # Actual lammps cell may be different from what is used to create R
        def fold(vec, pvec, i):
            p = pvec[i]
            x = vec[i] + 0.5 * p
            n = (np.mod(x, p) - x) / p
            return [float(self.f2qdec(a)) for a in (vec + n * pvec)]

        Apre[1, :] = fold(Apre[1, :], Apre[0, :], 0)
        Apre[2, :] = fold(Apre[2, :], Apre[1, :], 1)
        Apre[2, :] = fold(Apre[2, :], Apre[0, :], 0)

        self.A = Apre
        self.Ainv = np.linalg.inv(self.A)

        if self.is_skewed() and \
                (not (pbc[0] and pbc[1] and pbc[2])):
            raise RuntimeError('Skewed lammps cells MUST have '
                               'PBC == True in all directions!')

    def f2qdec(self, f):
        return dec.Decimal(repr(f)).quantize(self.car_prec, dec.ROUND_DOWN)

    def f2qs(self, f):
        return str(self.f2qdec(f))

    def f2s(self, f):
        return str(dec.Decimal(repr(f)).quantize(self.car_prec, dec.ROUND_HALF_EVEN))

    def dir2car(self, v):
        "Direct to cartesian coordinates"
        return np.dot(v, self.A)

    def car2dir(self, v):
        "Cartesian to direct coordinates"
        return np.dot(v, self.Ainv)

    def fold_to_str(self, v):
        "Fold a position into the lammps cell (semi open), return a tuple of str"
        # Two-stage fold, first into box, then into semi-open interval
        # (within the given precission).
        d = [x % (1 - self.dir_prec) for x in
             map(dec.Decimal, map(repr, np.mod(self.car2dir(v) + self.eps, 1.0)))]
        return tuple([self.f2qs(x) for x in
                      self.dir2car(list(map(float, d)))])

    def get_lammps_prism(self):
        A = self.A
        return (A[0, 0], A[1, 1], A[2, 2], A[1, 0], A[2, 0], A[2, 1])

    def get_lammps_prism_str(self):
        "Return a tuple of strings"
        p = self.get_lammps_prism()
        return tuple([self.f2s(x) for x in p])

    def pos_to_lammps_str(self, position):
        "Rotate an ase-cell position to the lammps cell orientation, return tuple of strs"
        return tuple([self.f2s(x) for x in np.dot(position, self.R)])

    def pos_to_lammps_fold_str(self, position):
        "Rotate and fold an ase-cell position into the lammps cell, return tuple of strs"
        return self.fold_to_str(np.dot(position, self.R))

    def is_skewed(self):
        acc = self.acc
        prism = self.get_lammps_prism()
        axy, axz, ayz = [np.abs(x) for x in prism[3:]]
        return (axy >= acc) or (axz >= acc) or (ayz >= acc)


def write_lammps_data(fileobj, atoms, specorder=None, force_skew=False,
                      prismobj=None, velocities=False):
    """Method which writes atomic structure data to a LAMMPS data file."""
    if isinstance(fileobj, str):
        f = paropen(fileobj, 'wb')
        close_file = True
    else:
        # Presume fileobj acts like a fileobj
        f = fileobj
        close_file = False

    if isinstance(atoms, list):
        if len(atoms) > 1:
            raise ValueError('Can only write one configuration to a lammps data file!')
        atoms = atoms[0]

    f.write('{0} (written by ASE) \n\n'.format(f.name).encode('utf-8'))

    symbols = atoms.get_chemical_symbols()
    n_atoms = len(symbols)
    f.write('{0} \t atoms \n'.format(n_atoms).encode('utf-8'))

    if specorder is None:
        # This way it is assured that LAMMPS atom types are always
        # assigned predictively according to the alphabetic order
        species = sorted(set(symbols))
    else:
        # To index elements in the LAMMPS data file
        # (indices must correspond to order in the potential file)
        species = specorder
    n_atom_types = len(species)
    f.write('{0}  atom types\n'.format(n_atom_types).encode('utf-8'))

    if prismobj is None:
        p = prism(atoms.get_cell())
    else:
        p = prismobj
    xhi, yhi, zhi, xy, xz, yz = p.get_lammps_prism_str()

    f.write('0.0 {0}  xlo xhi\n'.format(xhi).encode('utf-8'))
    f.write('0.0 {0}  ylo yhi\n'.format(yhi).encode('utf-8'))
    f.write('0.0 {0}  zlo zhi\n'.format(zhi).encode('utf-8'))

    if force_skew or p.is_skewed():
        f.write('{0} {1} {2}  xy xz yz\n'.format(xy, xz, yz).encode('utf-8'))
    f.write('\n\n'.encode('utf-8'))

    f.write('Atoms \n\n'.encode('utf-8'))
    for i, r in enumerate(map(p.pos_to_lammps_str,
                              atoms.get_positions())):
        s = species.index(symbols[i]) + 1
        f.write('{0:>6} {1:>3} {2} {3} {4}\n'.format(
            *(i + 1, s) + tuple(r)).encode('utf-8'))

    if velocities and atoms.get_velocities() is not None:
        f.write('\n\nVelocities \n\n'.encode('utf-8'))
        for i, v in enumerate(atoms.get_velocities()):
            f.write('{0:>6} {1} {2} {3}\n'.format(
                *(i + 1,) + tuple(v)).encode('utf-8'))

    f.flush()
    if close_file:
        f.close()


if __name__ == '__main__':
    pair_style = 'eam'
    Pd_eam_file = 'Pd_u3.eam'
    pair_coeff = ['* * ' + Pd_eam_file]
    parameters = {'pair_style': pair_style, 'pair_coeff': pair_coeff}
    files = [Pd_eam_file]
    calc = LAMMPS(parameters=parameters, files=files)
    a0 = 3.93
    b0 = a0 / 2.0
    if True:
        bulk = Atoms(['Pd'] * 4,
                     positions=[(0, 0, 0), (b0, b0, 0), (b0, 0, b0), (0, b0, b0)],
                     cell=[a0] * 3,
                     pbc=True)
        # test get_forces
        print('forces for a = {0}'.format(a0))
        print(calc.get_forces(bulk))
        # single points for various lattice constants
        bulk.set_calculator(calc)
        for n in range(-5, 5, 1):
            a = a0 * (1 + n / 100.0)
            bulk.set_cell([a] * 3)
            print('a : {0} , total energy : {1}'.format(a, bulk.get_potential_energy()))

    calc.clean()

