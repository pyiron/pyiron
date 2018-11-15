import numpy as np
import os
from subprocess import Popen, PIPE

from pyiron.vasp.outcar import Outcar
from pyiron.vasp.vasp import Vasp
from pyiron.vasp.structure import vasp_sorter
from pyiron.vasp.vasp import GenericOutput as GenericOutputBase
from pyiron.vasp.vasp import DFTOutput as DFTOutputBase
from pyiron.vasp.vasp import Output as OutputBase
from pyiron.interactive.generic import GenericInteractive


class VaspInt(GenericInteractive, Vasp):
    def __init__(self, project, job_name):
        super(VaspInt, self).__init__(project, job_name)
        self._interactive_write_input_files = True
        self._interactive_vasprun = None
        self.interactive_cache = {'cells': [],
                                  'energy_pot': [],
                                  'energy_tot': [],
                                  'forces': [],
                                  'positions': [],
                                  'indices': [],
                                  'steps': [],
                                  'computation_time': [],
                                  'volume': []}

    @property
    def interactive_enforce_structure_reset(self):
        return self._interactive_enforce_structure_reset

    @interactive_enforce_structure_reset.setter
    def interactive_enforce_structure_reset(self, reset):
        raise NotImplementedError('interactive_enforce_structure_reset() is not implemented!')

    def interactive_close(self):
        if self.interactive_is_activated():
            with open(os.path.join(self.working_directory, 'STOPCAR'), 'w') as stopcar:
                stopcar.write('LABORT = .TRUE.')  # stopcar.write('LSTOP = .TRUE.')
            try:
                self.run_if_interactive()
                self.run_if_interactive()
                for atom in self.current_structure.scaled_positions:
                    text = ' '.join(map('{:19.16f}'.format, atom))
                    self._interactive_library.stdin.write(text + '\n')
            except BrokenPipeError:
                self._logger.warn('VASP calculation exited before interactive_close() - already converged?')
            for key in self.interactive_cache.keys():
                if isinstance(self.interactive_cache[key], list):
                    self.interactive_cache[key] = self.interactive_cache[key][:-2]
            super(VaspInt, self).interactive_close()
            self.status.collect = True
            self._output_parser = Output()
            if self['vasprun.xml'] is not None:
                self.run()

    def interactive_energy_tot_getter(self):
        return self.interactive_energy_pot_getter()

    def interactive_energy_pot_getter(self):
        if self._interactive_vasprun is not None:
            file_name = os.path.join(self.working_directory, 'OUTCAR')
            return self._interactive_vasprun.get_energy_sigma_0(filename=file_name)[-1]
        else:
            return None

    def validate_ready_to_run(self):
        if self.server.run_mode.interactive and 'EDIFFG' in self.input.incar._dataset["Parameter"]:
            raise ValueError('If EDIFFG is defined VASP interrupts the interactive run_mode.')
        super(VaspInt, self).validate_ready_to_run()

    def interactive_forces_getter(self):
        if self._interactive_vasprun is not None:
            file_name = os.path.join(self.working_directory, 'OUTCAR')
            forces = self._interactive_vasprun.get_forces(filename=file_name)[-1]
            forces[vasp_sorter(self.structure)] = forces.copy()
            return forces
        else:
            return None

    def interactive_open(self):
        if self.executable.executable_path == '':
            self.status.aborted = True
            raise ValueError('No executable set!')
        if self.executable.mpi:
            self._interactive_library = Popen([self.executable.executable_path,
                                               str(self.server.cores)],
                                              stdout=PIPE,
                                              stdin=PIPE,
                                              stderr=PIPE,
                                              cwd=self.working_directory,
                                              universal_newlines=True)
        else:
            self._interactive_library = Popen(self.executable.executable_path,
                                              stdout=PIPE,
                                              stdin=PIPE,
                                              stderr=PIPE,
                                              cwd=self.working_directory,
                                              universal_newlines=True)

    def calc_minimize(self, e_tol=1e-8, f_tol=1e-8, max_iter=1000, pressure=None, n_print=1):
        raise NotImplementedError('calc_minimize() is not implemented for the interactive mode use calc_static()!')

    def calc_md(self, temperature=None, pressure=None, n_ionic_steps=1000, time_step=None, n_print=100, delta_temp=1.0,
                delta_press=None, seed=None, tloop=None, rescale_velocity=True):
        raise NotImplementedError('calc_md() is not implemented for the interactive mode use calc_static()!')

    def run_if_interactive(self):
        initial_run = not self.interactive_is_activated()
        super(VaspInt, self).run_if_interactive()
        if not initial_run:
            atom_numbers = self.current_structure.get_number_species_atoms()
            for species in atom_numbers.keys():
                indices = self.current_structure.select_index(species)
                for i in indices:
                    text = ' '.join(map('{:19.16f}'.format, self.current_structure.scaled_positions[i]))
                    self._logger.debug('Vasp library: ' + text)
                    self._interactive_library.stdin.write(text + '\n')
            self._interactive_library.stdin.flush()
        self._interactive_check_output()
        self._interactive_vasprun = Outcar()
        self.interactive_collect()

    def interactive_positions_setter(self, positions):
        pass

    def _check_incar_parameter(self, parameter, value):
        if parameter not in self.input.incar._dataset['Parameter']:
            self.input.incar[parameter] = value

    def _interactive_check_output(self):
        while self._interactive_library.poll() is None:
            text = self._interactive_library.stdout.readline()
            if "POSITIONS: reading from stdin" in text:
                return

    def _run_if_created(self, que_wait_for=None):
        if self.server.run_mode.interactive:
            self._check_incar_parameter(parameter='INTERACTIVE', value=True)
            self._check_incar_parameter(parameter='IBRION', value=-1)
            self._check_incar_parameter(parameter='POTIM', value=0.0)
            self._check_incar_parameter(parameter='NSW', value=1000)
            self._check_incar_parameter(parameter='ISYM', value=0)
        super(VaspInt, self)._run_if_created(que_wait_for=que_wait_for)


class Output(OutputBase):
    """
    Handles the output from a VASP simulation.

    Attributes:
        electronic_structure: Gives the electronic structure of the system
        electrostatic_potential: Gives the electrostatic/local potential of the system
        charge_density: Gives the charge density of the system
    """

    def __init__(self):
        super(Output, self).__init__()
        self.generic_output = GenericOutput()
        self.dft_output = DFTOutput()


class GenericOutput(GenericOutputBase):
    """

    This class stores the generic output like different structures, energies and forces from a simulation in a highly
    generic format. Usually the user does not have to access this class.

    Attributes:
        log_dict (dict): A dictionary of all tags and values of generic data (positions, forces, etc)
    """

    def __init__(self):
        super(GenericOutput, self).__init__()

    def to_hdf(self, hdf):
        """
        Save the object in a HDF5 file

        Args:
            hdf (pyiron_base.objects.generic.hdfio.ProjectHDFio): HDF path to which the object is to be saved

        """
        with hdf.open("generic") as hdf_go:
            # hdf_go["description"] = self.description
            for key, val in self.log_dict.items():
                if isinstance(val, list) or isinstance(val, np.ndarray):
                    hdf_go[key] = val[:-1]
                else:
                    hdf_go[key] = val
            with hdf_go.open("dft") as hdf_dft:
                for key, val in self.dft_log_dict.items():
                    if isinstance(val, list) or isinstance(val, np.ndarray):
                        hdf_dft[key] = val[:-1]
                    else:
                        hdf_dft[key] = val
                if self.bands.eigenvalue_matrix is not None:
                    self.bands.to_hdf_new(hdf_dft, "bands")


class DFTOutput(DFTOutputBase):
    """
    This class stores the DFT specific output

    Attributes:
        log_dict (dict): A dictionary of all tags and values of DFT data
    """

    def __init__(self):
        super(DFTOutput, self).__init__()

    def to_hdf(self, hdf):
        """
        Save the object in a HDF5 file

        Args:
            hdf (pyiron_base.objects.generic.hdfio.ProjectHDFio): HDF path to which the object is to be saved

        """
        with hdf.open("dft") as hdf_dft:
            # hdf_go["description"] = self.description
            for key, val in self.log_dict.items():
                if isinstance(val, list) or isinstance(val, np.ndarray):
                    hdf_dft[key] = val[:-1]
                else:
                    hdf_dft[key] = val
