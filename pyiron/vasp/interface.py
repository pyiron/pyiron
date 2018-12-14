from __future__ import print_function
import os
import posixpath
import numpy as np
import warnings
from shutil import copyfile
from subprocess import Popen, PIPE, check_output
from pyiron.atomistics.structure.atoms import Atoms
from pyiron.atomistics.md_analysis.trajectory_analysis import unwrap_coordinates
from pyiron.vasp.outcar import Outcar
from pyiron.vasp.procar import Procar
from pyiron.vasp.structure import read_atoms, vasp_sorter
from pyiron.vasp.vasprun import Vasprun as Vr
from pyiron.vasp.vasprun import VasprunError
from pyiron.vasp.volumetric_data import VaspVolumetricData
from pyiron.dft.waves.electronic import ElectronicStructure
from pyiron.base.job.interface import FileInterface
from pyiron.atomistics.job.interface import AtomisticInteractiveInterface


class VaspInterface(FileInterface):
    def __init__(self):
        self.output_parser = Output()

    def write_input(self, job):
        """
        Call routines that generate the INCAR, POTCAR, KPOINTS and POSCAR input files
        """
        if job.input.incar['SYSTEM'] == 'pyiron_jobname':
            job.input.incar['SYSTEM'] = job.job_name
        job.write_magmoms()
        job.set_coulomb_interactions()
        if "CONTCAR" in job.restart_file_dict.keys():
            if job.restart_file_dict["CONTCAR"] == "POSCAR":
                if job.server.run_mode.modal:
                    warnings.warn(
                        "The POSCAR file will be overwritten by the CONTCAR file specified in restart_file_list.")
                else:
                    job.logger.info(
                        "The POSCAR file will be overwritten by the CONTCAR file specified in restart_file_list.")
        job.input.write(structure=job.structure, directory=job.working_directory)

    # define routines that collect all output files
    def collect_output(self, job):
        """
        Collects the outputs and stores them to the hdf file
        """
        if job.structure is None or len(job.structure) == 0:
            try:
                job.structure = self.get_final_structure_from_file(job=job, filename="CONTCAR")
            except IOError:
                job.structure = self.get_final_structure_from_file(job=job, filename="POSCAR")
        self.output_parser.structure = job.structure.copy()
        try:
            self.output_parser.collect(directory=job.working_directory)
        except VaspCollectError:
            job.status.aborted = True
            return
        self.output_parser.to_hdf(job._hdf5)

    def collect_logfiles(self, job):
        """
        Collect errors and warnings.
        """
        self.collect_errors(job=job)
        self.collect_warnings(job=job)

    def collect_warnings(self, job):
        """
        Collects warnings from the VASP run
        """
        # TODO: implement for VASP
        job._logger.info("collect_warnings() is not yet implemented for VASP")

    def collect_errors(self, job):
        """
        Collects errors from the VASP run
        """
        # TODO: implement for vasp
        job._logger.info("collect_errors() is not yet implemented for VASP")

    def get_final_structure_from_file(self, job, filename="CONTCAR"):
        """
        Get the final structure of the simulation usually from the CONTCAR file

        Args:
            filename (str): Path to the CONTCAR file in VASP

        Returns:
            pyiron.atomistics.structure.atoms.Atoms: The final structure
        """
        filename = posixpath.join(job.working_directory, filename)
        input_structure = job.structure.copy()
        try:
            output_structure = read_atoms(filename=filename, species_list=input_structure.get_parent_elements())
        except (IndexError, ValueError, IOError):
            job._logger.warning("Unable to read output structure")
            return
        input_structure.cell = output_structure.cell.copy()
        input_structure.positions[job.sorted_indices] = output_structure.positions
        return input_structure

    def cleanup(self, job, files_to_remove=("WAVECAR", "CHGCAR", "CHG", "vasprun.xml")):
        """
        Removes excess files (by default: WAVECAR, CHGCAR, CHG)
        """
        list_files = job.list_files()
        for file in list_files:
            if file in files_to_remove:
                abs_file_path = os.path.join(job.working_directory, file)
                os.remove(abs_file_path)

    def from_directory(self, job, directory):
        """
        The Vasp instance is created by parsing the input and output from the specified directory

        Args:
            directory (str): Path to the directory
        """
        if not job.status.finished:
            # _ = s.top_path(directory)
            files = os.listdir(directory)
            vp_new = Vr()
            if "OUTCAR.gz" in files and "OUTCAR" not in files:
                _ = check_output(['gzip', '-d', 'OUTCAR.gz'], cwd=directory, shell=False, universal_newlines=True)
                files = os.listdir(directory)
            if "vasprun.xml.bz2" in files:
                _ = check_output(['bzip2', '-d', 'vasprun.xml.bz2'], cwd=directory, shell=False, universal_newlines=True)
                files = os.listdir(directory)
            if "vasprun.xml.gz" in files:
                _ = check_output(['gzip', '-d', 'vasprun.xml.gz'], cwd=directory, shell=False, universal_newlines=True)
                files = os.listdir(directory)
            try:
                if not ("OUTCAR" in files or "vasprun.xml" in files):
                    raise IOError("This file isn't present")
                    # raise AssertionError("OUTCAR/vasprun.xml should be present in order to import from directory")
                if "vasprun.xml" in files:
                    vp_new.from_file(filename=posixpath.join(directory, "vasprun.xml"))
                    job.structure = vp_new.get_initial_structure()
            except (IOError, VasprunError):  # except AssertionError:
                pass
                # raise AssertionError("OUTCAR/vasprun.xml should be present in order to import from directory")
            if "INCAR" in files:
                try:
                    job.input.incar.read_input(posixpath.join(directory, "INCAR"), ignore_trigger="!")
                except (IndexError, TypeError, ValueError):
                    pass
            if "KPOINTS" in files:
                try:
                    job.input.kpoints.read_input(posixpath.join(directory, "KPOINTS"), ignore_trigger="!")
                except (IndexError, TypeError, ValueError):
                    pass
            if "POSCAR" in files and "POTCAR" in files:
                structure = read_atoms(posixpath.join(directory, "POSCAR"), species_from_potcar=True)
            else:
                structure = vp_new.get_initial_structure()
            job.structure = structure
            job._write_chemical_formular_to_database()
            job._import_directory = directory
            job.status.collect = True
            job.to_hdf()
            job.collect_output()
            job.status.finished = True
        else:
            return

    @staticmethod
    def stop_calculation(job, next_electronic_step=False):
        """
        Call to stop the VASP calculation

        Args:
            next_electronic_step (bool): True if the next electronic step should be calculated

        """
        filename = os.path.join(job.working_directory, 'STOPCAR')
        with open(filename, 'w') as f:
            if not next_electronic_step:
                f.write('LSTOP = .TRUE.\n')
            else:
                f.write('LABORT =.TRUE.\n')

    def reset_output(self):
        """
        Resets the output instance
        """
        self.output_parser = Output()

    def copy_chgcar(self, job, old_vasp_job):
        """
        Copy CHGCAR from previous VASP calcualtion to the new VASP job.
        (Sets ICHARG = 1)

        Args:
            old_vasp_job (pyiron.vasp.vasp.Vasp): Finished Vasp job instance

        """
        self.copy_file(job=job, old_vasp_job=old_vasp_job)
        job.input.incar["ICHARG"] = 1

    def copy_wavecar(self, job, old_vasp_job):
        """
        Copy WAVECAR from previous VASP calculation to the new VASP job.
        (Sets ICHARG = 1)

        Args:
            (pyiron.vasp.vasp.Vasp): Finished Vasp job instance

        """
        self.copy_file(job=job, old_vasp_job=old_vasp_job, filename="WAVECAR")
        job.input.incar["ISTART"] = 1

    def copy_file(self, job, old_vasp_job, filename="CHGCAR"):
        """
        Copy a file from a previous vasp job

        Args:
            old_vasp_job (pyiron.vasp.vasp.Vasp): Finished Vasp job instance
            filename (str): Destination to copy the file

        """
        old_path = os.path.join(old_vasp_job.working_directory, filename)
        new_path = os.path.join(job.working_directory, filename)
        if not os.path.isdir(job.working_directory):
            os.makedirs(job.working_directory)
        copyfile(old_path, new_path)


class InteractiveVaspInterface(VaspInterface, AtomisticInteractiveInterface):
    def __init__(self):
        super(InteractiveVaspInterface, self).__init__()
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

    def interactive_close(self, job):
        if self.interactive_is_activated():
            with open(os.path.join(job.working_directory, 'STOPCAR'), 'w') as stopcar:
                stopcar.write('LABORT = .TRUE.')  # stopcar.write('LSTOP = .TRUE.')
            try:
                self.run_if_interactive(job=job)
                self.run_if_interactive(job=job)
                for atom in job.current_structure.scaled_positions:
                    text = ' '.join(map('{:19.16f}'.format, atom))
                    self._interactive_library.stdin.write(text + '\n')
            except (BrokenPipeError, IOError):
                job._logger.warn('VASP calculation exited before interactive_close() - already converged?')
            for key in self.interactive_cache.keys():
                if isinstance(self.interactive_cache[key], list):
                    self.interactive_cache[key] = self.interactive_cache[key][:-2]
            super(InteractiveVaspInterface, self).interactive_close(job=job)
            job.status.collect = True
            self._output_parser = Output()
            if job['vasprun.xml'] is not None:
                job.run()

    def interactive_energy_tot_getter(self, job):
        return self.interactive_energy_pot_getter(job=job)

    def interactive_energy_pot_getter(self, job):
        if self._interactive_vasprun is not None:
            file_name = os.path.join(job.working_directory, 'OUTCAR')
            return self._interactive_vasprun.get_energy_sigma_0(filename=file_name)[-1]
        else:
            return None

    def interactive_forces_getter(self, job):
        if self._interactive_vasprun is not None:
            file_name = os.path.join(job.working_directory, 'OUTCAR')
            forces = self._interactive_vasprun.get_forces(filename=file_name)[-1]
            forces[vasp_sorter(job.structure)] = forces.copy()
            return forces
        else:
            return None

    def interactive_open(self, job):
        if job.executable.executable_path == '':
            job.status.aborted = True
            raise ValueError('No executable set!')
        if job.executable.mpi:
            self._interactive_library = Popen([job.executable.executable_path,
                                               str(job.server.cores)],
                                              stdout=PIPE,
                                              stdin=PIPE,
                                              stderr=PIPE,
                                              cwd=job.working_directory,
                                              universal_newlines=True)
        else:
            self._interactive_library = Popen(job.executable.executable_path,
                                              stdout=PIPE,
                                              stdin=PIPE,
                                              stderr=PIPE,
                                              cwd=job.working_directory,
                                              universal_newlines=True)

    def run_if_interactive_non_modal(self, job):
        initial_run = not self.interactive_is_activated()
        super(InteractiveVaspInterface, self).run_if_interactive()
        if not initial_run:
            atom_numbers = job.current_structure.get_number_species_atoms()
            for species in atom_numbers.keys():
                indices = job.current_structure.select_index(species)
                for i in indices:
                    text = ' '.join(map('{:19.16f}'.format, job.current_structure.scaled_positions[i]))
                    job._logger.debug('Vasp library: ' + text)
                    self._interactive_library.stdin.write(text + '\n')
            self._interactive_library.stdin.flush()
        self._interactive_fetch_completed = False

    def run_if_interactive(self, job):
        self.run_if_interactive_non_modal()
        self._interactive_check_output()
        self._interactive_vasprun = Outcar()
        self.interactive_collect()

    def interactive_fetch(self, job):
        if self._interactive_fetch_completed and job.server.run_mode.interactive_non_modal:
            print('First run and then fetch')
        else:
            self._interactive_check_output()
            self._interactive_vasprun = Outcar()
            super(InteractiveVaspInterface, self).interactive_collect()
            job._logger.debug('interactive run - done')

    def interactive_positions_setter(self, job, positions):
        pass

    def _check_incar_parameter(self, job, parameter, value):
        if parameter not in job.input.incar._dataset['Parameter']:
            job.input.incar[parameter] = value

    def _interactive_check_output(self):
        while self._interactive_library.poll() is None:
            text = self._interactive_library.stdout.readline()
            if "POSITIONS: reading from stdin" in text:
                return


class Output(object):
    """
    Handles the output from a VASP simulation.

    Attributes:
        electronic_structure: Gives the electronic structure of the system
        electrostatic_potential: Gives the electrostatic/local potential of the system
        charge_density: Gives the charge density of the system
    """

    def __init__(self):
        self._structure = None
        self.outcar = Outcar()
        self.generic_output = GenericOutput()
        self.description = "This contains all the output static from this particular vasp run"
        self.charge_density = VaspVolumetricData()
        self.electrostatic_potential = VaspVolumetricData()
        self.procar = Procar()
        self.electronic_structure = ElectronicStructure()
        self.vp_new = Vr()

    @property
    def structure(self):
        """
        Getter for the output structure
        """
        return self._structure

    @structure.setter
    def structure(self, atoms):
        """
        Setter for the output structure
        """
        self._structure = atoms

    def collect(self, directory=os.getcwd()):
        """
        Collects output from the working directory

        Args:
            directory (str): Path to the directory
        """
        sorted_indices = vasp_sorter(self.structure)
        files_present = os.listdir(directory)
        log_dict = dict()
        vasprun_working, outcar_working = False, False
        if not ("OUTCAR" in files_present or "vasprun.xml" in files_present):
            raise IOError("Either the OUTCAR or vasprun.xml files need to be present")
        if "OUTCAR" in files_present:
            self.outcar.from_file(filename=posixpath.join(directory, "OUTCAR"))
            outcar_working = True
        if "vasprun.xml" in files_present:
            try:
                self.vp_new.from_file(filename=posixpath.join(directory, "vasprun.xml"))
            except VasprunError:
                pass
            else:
                vasprun_working = True

        if outcar_working:
            log_dict["temperature"] = self.outcar.parse_dict["temperatures"]
            log_dict["pressures"] = self.outcar.parse_dict["pressures"]
            self.generic_output.dft_log_dict["n_elect"] = self.outcar.parse_dict["n_elect"]
            if len(self.outcar.parse_dict["magnetization"]) > 0:
                magnetization = np.array(self.outcar.parse_dict["magnetization"]).copy()
                final_magmoms = np.array(self.outcar.parse_dict["final_magmoms"]).copy()
                # magnetization[sorted_indices] = magnetization.copy()
                if len(final_magmoms) != 0:
                    if len(final_magmoms.shape) == 3:
                        final_magmoms[:, sorted_indices, :] = final_magmoms.copy()
                    else:
                        final_magmoms[:, sorted_indices] = final_magmoms.copy()
                self.generic_output.dft_log_dict["magnetization"] = magnetization.tolist()
                self.generic_output.dft_log_dict["final_magmoms"] = final_magmoms.tolist()

        if vasprun_working:
            log_dict["forces"] = self.vp_new.vasprun_dict["forces"]
            log_dict["cells"] = self.vp_new.vasprun_dict["cells"]
            log_dict["volume"] = [np.linalg.det(cell) for cell in self.vp_new.vasprun_dict["cells"]]
            # log_dict["total_energies"] = self.vp_new.vasprun_dict["total_energies"]
            log_dict["energy_tot"] = self.vp_new.vasprun_dict["total_energies"]
            if "kinetic_energies" in self.vp_new.vasprun_dict.keys():
                log_dict["energy_pot"] = log_dict["energy_tot"] - self.vp_new.vasprun_dict["kinetic_energies"]
            else:
                log_dict["energy_pot"] = log_dict["energy_tot"]
            log_dict["steps"] = np.arange(len(log_dict["energy_tot"]))
            log_dict["positions"] = self.vp_new.vasprun_dict["positions"]
            log_dict["forces"][:, sorted_indices] = log_dict["forces"].copy()
            log_dict["positions"][:, sorted_indices] = log_dict["positions"].copy()
            log_dict["positions_unwrapped"] = unwrap_coordinates(positions=log_dict["positions"], cell=None,
                                                                 is_relative=True)
            for i, pos in enumerate(log_dict["positions"]):
                log_dict["positions"][i] = np.dot(pos, log_dict["cells"][i])
                log_dict["positions_unwrapped"][i] = np.dot(log_dict["positions_unwrapped"][i].copy(),
                                                            log_dict["cells"][i])
            # log_dict["scf_energies"] = self.vp_new.vasprun_dict["scf_energies"]
            # log_dict["scf_dipole_moments"] = self.vp_new.vasprun_dict["scf_dipole_moments"]
            self.electronic_structure = self.vp_new.get_electronic_structure()
            if self.electronic_structure.grand_dos_matrix is not None:
                self.electronic_structure.grand_dos_matrix[:, :, :, sorted_indices, :] = \
                    self.electronic_structure.grand_dos_matrix[:, :, :, :, :].copy()
            if self.electronic_structure.resolved_densities is not None:
                self.electronic_structure.resolved_densities[:, sorted_indices, :, :] = \
                    self.electronic_structure.resolved_densities[:, :, :, :].copy()
            self.structure.positions = log_dict["positions"][-1]
            self.structure.cell = log_dict["cells"][-1]

        elif outcar_working:
            # log_dict = self.outcar.parse_dict.copy()
            if len(self.outcar.parse_dict["energies"]) == 0:
                raise VaspCollectError("Error in parsing OUTCAR")
            log_dict["energy_tot"] = self.outcar.parse_dict["energies"]
            log_dict["temperature"] = self.outcar.parse_dict["temperatures"]
            log_dict["pressures"] = self.outcar.parse_dict["pressures"]
            log_dict["forces"] = self.outcar.parse_dict["forces"]
            log_dict["positions"] = self.outcar.parse_dict["positions"]
            # log_dict["forces"][:, sorted_indices] = log_dict["forces"].copy()
            # log_dict["positions"][:, sorted_indices] = log_dict["positions"].copy()
            if len(log_dict["positions"].shape) != 3:
                raise VaspCollectError("Improper OUTCAR parsing")
            elif log_dict["positions"].shape[1] != len(sorted_indices):
                raise VaspCollectError("Improper OUTCAR parsing")
            if len(log_dict["forces"].shape) != 3:
                raise VaspCollectError("Improper OUTCAR parsing")
            elif log_dict["forces"].shape[1] != len(sorted_indices):
                raise VaspCollectError("Improper OUTCAR parsing")
            log_dict["time"] = self.outcar.parse_dict["time"]
            log_dict["steps"] = self.outcar.parse_dict["steps"]
            log_dict["cells"] = self.outcar.parse_dict["cells"]
            log_dict["volume"] = np.array([np.linalg.det(cell) for cell in self.outcar.parse_dict["cells"]])
            self.generic_output.dft_log_dict["scf_energy_free"] = self.outcar.parse_dict["scf_energies"]
            self.generic_output.dft_log_dict["scf_dipole_mom"] = self.outcar.parse_dict["scf_dipole_moments"]
            self.generic_output.dft_log_dict["n_elect"] = self.outcar.parse_dict["n_elect"]
            self.generic_output.dft_log_dict["energy_int"] = self.outcar.parse_dict["energies_int"]
            self.generic_output.dft_log_dict["energy_free"] = self.outcar.parse_dict["energies"]
            self.generic_output.dft_log_dict["energy_zero"] = self.outcar.parse_dict["energies_zero"]
            if "PROCAR" in files_present:
                try:
                    self.electronic_structure = self.procar.from_file(filename=posixpath.join(directory, "PROCAR"))
                    #  Even the atom resolved values have to be sorted from the vasp atoms order to the Atoms order
                    self.electronic_structure.grand_dos_matrix[:, :, :, sorted_indices, :] = \
                        self.electronic_structure.grand_dos_matrix[:, :, :, :, :].copy()
                    try:
                        self.electronic_structure.efermi = self.outcar.parse_dict["fermi_level"]
                    except KeyError:
                        self.electronic_structure.efermi = self.vp_new.vasprun_dict["efermi"]
                except ValueError:
                    pass

        # important that we "reverse sort" the atoms in the vasp format into the atoms in the atoms class
        self.generic_output.log_dict = log_dict
        if vasprun_working:
            # self.dft_output.log_dict["parameters"] = self.vp_new.vasprun_dict["parameters"]
            self.generic_output.dft_log_dict["scf_dipole_mom"] = self.vp_new.vasprun_dict["scf_dipole_moments"]
            if len(self.generic_output.dft_log_dict["scf_dipole_mom"][0]) > 0:
                total_dipole_moments = np.array([dip[-1] for dip in self.generic_output.dft_log_dict["scf_dipole_mom"]])
                self.generic_output.dft_log_dict["dipole_mom"] = total_dipole_moments
            self.generic_output.dft_log_dict["scf_energy_int"] = self.vp_new.vasprun_dict["scf_energies"]
            self.generic_output.dft_log_dict["scf_energy_free"] = self.vp_new.vasprun_dict["scf_fr_energies"]
            self.generic_output.dft_log_dict["scf_energy_zero"] = self.vp_new.vasprun_dict["scf_0_energies"]
            self.generic_output.dft_log_dict["energy_int"] = np.array([e_int[-1] for e_int in
                                                                      self.generic_output.dft_log_dict
                                                                      ["scf_energy_int"]])
            self.generic_output.dft_log_dict["energy_free"] = np.array([e_free[-1] for e_free in
                                                                       self.generic_output.dft_log_dict
                                                                       ["scf_energy_free"]])
            self.generic_output.dft_log_dict["energy_zero"] = np.array([e_zero[-1] for e_zero in
                                                                       self.generic_output.dft_log_dict
                                                                       ["scf_energy_zero"]])
            self.generic_output.dft_log_dict["n_elect"] = float(self.vp_new.vasprun_dict["parameters"]["electronic"]
                                                                ['NELECT'])
            if "kinetic_energies" in self.vp_new.vasprun_dict.keys():
                self.generic_output.dft_log_dict["scf_energy_kin"] = self.vp_new.vasprun_dict["kinetic_energies"]

        if "LOCPOT" in files_present:
            self.electrostatic_potential.from_file(filename=posixpath.join(directory, "LOCPOT"), normalize=False)
        if "CHGCAR" in files_present:
            self.charge_density.from_file(filename=posixpath.join(directory, "CHGCAR"), normalize=True)
        self.generic_output.bands = self.electronic_structure

    def to_hdf(self, hdf):
        """
        Save the object in a HDF5 file

        Args:
            hdf (pyiron.base.generic.hdfio.ProjectHDFio): HDF path to which the object is to be saved

        """
        with hdf.open("output") as hdf5_output:
            hdf5_output["description"] = self.description
            self.generic_output.to_hdf(hdf5_output)
            try:
                self.structure.to_hdf(hdf5_output)
            except AttributeError:
                pass

            # with hdf5_output.open("vasprun") as hvr:
            #  if self.vasprun.dict_vasprun is not None:
            #     for key, val in self.vasprun.dict_vasprun.items():
            #        hvr[key] = val

            if self.electrostatic_potential.total_data is not None:
                self.electrostatic_potential.to_hdf(hdf5_output, group_name="electrostatic_potential")

            if self.charge_density.total_data is not None:
                self.charge_density.to_hdf(hdf5_output, group_name="charge_density")

            if len(self.electronic_structure.kpoint_list) > 0:
                self.electronic_structure.to_hdf(hdf=hdf5_output, group_name="electronic_structure")

            if self.outcar.parse_dict:
                self.outcar.to_hdf_minimal(hdf=hdf5_output, group_name="outcar")

    def from_hdf(self, hdf):
        """
        Reads the attributes and reconstructs the object from a hdf file
        Args:
            hdf: The hdf5 instance
        """
        with hdf.open("output") as hdf5_output:
            # self.description = hdf5_output["description"]
            if self.structure is None:
                self.structure = Atoms()
            self.structure.from_hdf(hdf5_output)
            self.generic_output.from_hdf(hdf5_output)
            if "electrostatic_potential" in hdf5_output.list_groups():
                self.electrostatic_potential.from_hdf(hdf5_output, group_name="electrostatic_potential")
            if "charge_density" in hdf5_output.list_groups():
                self.charge_density.from_hdf(hdf5_output, group_name="charge_density")
            if "electronic_structure" in hdf5_output.list_groups():
                self.electronic_structure.from_hdf(hdf=hdf5_output)
            if "outcar" in hdf5_output.list_groups():
                self.outcar.from_hdf(hdf=hdf5_output, group_name="outcar")


class GenericOutput(object):
    """

    This class stores the generic output like different structures, energies and forces from a simulation in a highly
    generic format. Usually the user does not have to access this class.

    Attributes:
        log_dict (dict): A dictionary of all tags and values of generic data (positions, forces, etc)
    """

    def __init__(self):
        self.log_dict = dict()
        self.dft_log_dict = dict()
        self.description = "generic_output contains generic output static"
        self._bands = ElectronicStructure()

    @property
    def bands(self):
        return self._bands

    @bands.setter
    def bands(self, val):
        self._bands = val

    def to_hdf(self, hdf):
        """
        Save the object in a HDF5 file

        Args:
            hdf (pyiron.base.generic.hdfio.ProjectHDFio): HDF path to which the object is to be saved

        """
        with hdf.open("generic") as hdf_go:
            # hdf_go["description"] = self.description
            for key, val in self.log_dict.items():
                hdf_go[key] = val
            with hdf_go.open("dft") as hdf_dft:
                for key, val in self.dft_log_dict.items():
                    hdf_dft[key] = val
                if self.bands.eigenvalue_matrix is not None:
                    self.bands.to_hdf_new(hdf_dft, "bands")

    def from_hdf(self, hdf):
        """
        Reads the attributes and reconstructs the object from a hdf file
        Args:
            hdf: The hdf5 instance
        """
        with hdf.open("generic") as hdf_go:
            for node in hdf_go.list_nodes():
                if node == "description":
                    # self.description = hdf_go[node]
                    pass
                else:
                    self.log_dict[node] = hdf_go[node]
            if 'dft' in hdf_go.list_groups():
                with hdf_go.open("dft") as hdf_dft:
                    for node in hdf_dft.list_nodes():
                        self.dft_log_dict[node] = hdf_dft[node]
                    if 'bands' in hdf_dft.list_groups():
                        self.bands.from_hdf_new(hdf_dft, "bands")


class VaspCollectError(ValueError):
    pass
