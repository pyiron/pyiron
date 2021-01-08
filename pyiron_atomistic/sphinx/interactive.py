# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
import os
import scipy.constants
import subprocess
import warnings
import time
from pyiron_atomistic.sphinx.base import SphinxBase, Group
from pyiron_atomistic.atomistics.job.interactive import GenericInteractive, GenericInteractiveOutput
from pyiron_atomistic.vasp.potential import VaspPotentialSetter

BOHR_TO_ANGSTROM = (
    scipy.constants.physical_constants["Bohr radius"][0] / scipy.constants.angstrom
)
HARTREE_TO_EV = scipy.constants.physical_constants["Hartree energy in eV"][0]
HARTREE_OVER_BOHR_TO_EV_OVER_ANGSTROM = HARTREE_TO_EV / BOHR_TO_ANGSTROM

__author__ = "Osamu Waseda, Jan Janssen"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "development"
__date__ = "Sep 1, 2017"


class SphinxInteractive(SphinxBase, GenericInteractive):
    def __init__(self, project, job_name):
        super(SphinxInteractive, self).__init__(project, job_name)
        self.output = SphinxOutput(job=self)
        self._interactive_write_input_files = True
        self._interactive_library_read = None
        self._interactive_fetch_completed = True
        self.interactive_flush_frequency = 1

    @property
    def structure(self):
        return GenericInteractive.structure.fget(self)

    @structure.setter
    def structure(self, structure):
        GenericInteractive.structure.fset(self, structure)
        if structure is not None:
            self._potential = VaspPotentialSetter(
                element_lst=structure.get_species_symbols().tolist()
            )

    def get_structure(self, iteration_step=-1, wrap_atoms=True):
        return GenericInteractive.get_structure(
            self, iteration_step=iteration_step, wrap_atoms=wrap_atoms
        )

    def interactive_energy_tot_getter(self):
        return self.interactive_energy_pot_getter()

    def interactive_energy_pot_getter(self):
        self._interactive_pipe_write("get energy")
        return float(self._interactive_library_read.readline()) * HARTREE_TO_EV

    def interactive_forces_getter(self):
        self._interactive_pipe_write("get forces")
        ff = []
        for _ in range(len(self.structure)):
            line = self._interactive_library_read.readline().split()
            ff.append(
                [
                    float(line[i]) * HARTREE_OVER_BOHR_TO_EV_OVER_ANGSTROM
                    for i in range(3)
                ]
            )
        ff = np.array(ff)[self.id_spx_to_pyi]
        return ff

    def interactive_cells_getter(self):
        self._interactive_pipe_write("get cell")
        cc = []
        for _ in range(3):
            line = self._interactive_library_read.readline().split()
            cc.append([float(line[i]) * BOHR_TO_ANGSTROM for i in range(3)])
        return np.array(cc)

    def interactive_positions_getter(self):
        self._interactive_pipe_write("get structure")
        xx = []
        for _ in range(len(self.structure)):
            line = self._interactive_library_read.readline().split()
            xx.append([float(line[i]) * BOHR_TO_ANGSTROM for i in range(3)])
        xx = np.array(xx)[self.id_spx_to_pyi]
        return xx

    def interactive_positions_setter(self, positions):
        self._interactive_pipe_write("set structure")
        positions = positions[self.id_pyi_to_spx]
        positions = np.reshape(positions, 3 * len(self.structure)) / BOHR_TO_ANGSTROM
        self._interactive_pipe_write(positions.tolist())

    def interactive_spins_getter(self):
        self._logger.debug("get spins - start ...")
        self._interactive_pipe_write("get atomspin")
        mm = []
        for _ in range(len(self.structure)):
            line = self._interactive_library_read.readline().split()
            mm.append(float(line[0]))
        mm = np.array(mm)[self.id_spx_to_pyi]
        # self.interactive_cache['atom_spins'].append(mm)
        self._logger.debug("get spins - done.")
        return mm

    def interactive_spin_constraints_setter(self, spins):
        if self._generic_input["fix_spin_constraint"]:
            self._logger.debug("set spin constraints - start ...")
            self._spin_constraint_enabled = True
            self._interactive_pipe_write("set spinconstraint")
            spins = np.array(spins)[self.id_pyi_to_spx]
            self._spin_constraints = np.array(spins)
            self._interactive_pipe_write(spins.tolist())
            # self.interactive_cache['atom_spin_constraints'].append(spins)
            self._logger.debug("set spin constraints - done.")
        else:
            warnings.warn("Spin constraint not set -> set fix_spin_constraint = True")

    def interactive_spin_constraints_getter(self):
        return self._spin_constraints
        # return self.interactive_cache['atom_spin_constraints'][-1]

    def interactive_magnetic_forces_getter(self):
        if self._generic_input["fix_spin_constraint"]:
            self._interactive_pipe_write("get nu")
            nn = []
            for _ in range(len(self.structure)):
                line = self._interactive_library_read.readline().split()
                nn.append(HARTREE_TO_EV * float(line[0]))
            nn = np.array(nn)[self.id_spx_to_pyi]
            return nn
        else:
            return None

    def interactive_initialize_interface(self):
        self.server.threads = self.input["THREADS"]
        if self.executable.executable_path == "":
            self.status.aborted = True
            raise ValueError("No executable set!")
        if self.server.cores == 1 or not self.executable.mpi:
            out = subprocess.Popen(
                str(self.executable),
                cwd=self.project_hdf5.working_directory,
                shell=True,
                stderr=subprocess.STDOUT,
                universal_newlines=True,
            )
        else:
            out = subprocess.Popen(
                [
                    self.executable.executable_path,
                    str(self.server.cores),
                    str(self.server.threads),
                ],
                cwd=self.project_hdf5.working_directory,
                shell=False,
                stderr=subprocess.STDOUT,
                universal_newlines=True,
            )
        while not self._interactive_pipes_initialized:
            time.sleep(1)
        self._logger.debug("open interactive interface!")
        self._interactive_library = open(
            os.path.join(self.working_directory, "sxctrl"), "w"
        )
        self._interactive_library_read = open(
            os.path.join(self.working_directory, "sxres"), "r"
        )
        self._logger.debug("interactive interface is opened!")
        if (
            np.all(self.structure.get_initial_magnetic_moments() == None)
            and "atom_spins" in self.interactive_cache.keys()
        ):
            del self.interactive_cache["atom_spins"]
        if self._generic_input["fix_spin_constraint"]:
            self.interactive_spin_constraints_setter(
                self._structure_current.get_initial_magnetic_moments()
            )
        else:
            if "magnetic_forces" in self.interactive_cache.keys():
                del self.interactive_cache["magnetic_forces"]
            if "atom_spin_constraints" in self.interactive_cache.keys():
                del self.interactive_cache["atom_spin_constraints"]
        if len(self.restart_file_list) > 0:
            self._logger.debug("restarting; interactive run - start ...")
            self._interactive_pipe_write("run restart")
            self.interactive_fetch()

    def _output_interactive_to_generic(self):
        with self.project_hdf5.open("output") as h5:
            if "interactive" in h5.list_groups():
                for key in ["positions", "cells", "indices", "cells", "forces"]:
                    h5["generic/" + key] = h5["interactive/" + key]
                with self.project_hdf5.open("input") as hdf5_input:
                    with hdf5_input.open("generic") as hdf5_generic:
                        if "dft" not in hdf5_generic.list_groups():
                            hdf5_generic.create_group("dft")
                        with hdf5_generic.open("dft") as hdf5_dft:
                            if (
                                "atom_spin_constraints"
                                in h5["interactive"].list_nodes()
                            ):
                                hdf5_dft["atom_spin_constraints"] = h5[
                                    "interactive/atom_spin_constraints"
                                ]

    def collect_output(self, force_update=False):
        super(SphinxInteractive, self).collect_output(force_update=force_update)
        self._output_interactive_to_generic()

    def interactive_close(self):
        if self.interactive_is_activated():
            self._interactive_pipe_write("end")
            self._interactive_library.close()
            self._interactive_library_read.close()
            self.status.collect = True
            if self["energy.dat"] is not None:
                self.run()
            self._output_interactive_to_generic()
            super(SphinxInteractive, self).interactive_close()

    def calc_minimize(
        self,
        electronic_steps=None,
        ionic_steps=None,
        max_iter=None,
        pressure=None,
        algorithm=None,
        retain_charge_density=False,
        retain_electrostatic_potential=False,
        ionic_energy_tolerance=0.0,
        ionic_force_tolerance=1.0e-2,
        ionic_energy=None,
        ionic_forces=None,
        volume_only=False,
    ):
        if (
            self.server.run_mode.interactive
            or self.server.run_mode.interactive_non_modal
        ):
            raise NotImplementedError(
                "calc_minimize() is not implemented for the interactive mode use calc_static()!"
            )
        else:
            super(SphinxInteractive, self).calc_minimize(
                electronic_steps=electronic_steps,
                ionic_steps=ionic_steps,
                max_iter=max_iter,
                pressure=pressure,
                algorithm=algorithm,
                retain_charge_density=retain_charge_density,
                retain_electrostatic_potential=retain_electrostatic_potential,
                ionic_energy_tolerance=ionic_energy_tolerance,
                ionic_force_tolerance=ionic_force_tolerance,
                volume_only=volume_only,
            )

    def run_if_interactive(self):
        super(SphinxInteractive, self).run_if_interactive()
        self._logger.debug("interactive run - start ...")
        self._interactive_pipe_write("run electronicminimization")
        self.interactive_fetch()

    def run_if_interactive_non_modal(self):
        if not self._interactive_fetch_completed:
            print("Warning: interactive_fetch being effectuated")
            self.interactive_fetch()
        super(SphinxInteractive, self).run_if_interactive()
        self._logger.debug("interactive run - start ...")
        self._interactive_pipe_write("run electronicminimization")
        self._interactive_fetch_completed = False

    def interactive_fetch(self):
        if (
            self._interactive_fetch_completed
            and self.server.run_mode.interactive_non_modal
        ):
            print("First run and then fetch")
        else:
            self.interactive_collect()
            self._logger.debug("interactive run - done")

    @property
    def _interactive_pipes_initialized(self):
        return os.path.exists(
            os.path.join(self.working_directory, "sxctrl")
        ) and os.path.exists(os.path.join(self.working_directory, "sxres"))

    def _interactive_pipe_write(self, line):
        if isinstance(line, str) or isinstance(line, int) or isinstance(line, float):
            self._interactive_library.write(str(line) + "\n")
            self._interactive_library.flush()
        elif isinstance(line, list):
            for subline in line:
                self._interactive_pipe_write(subline)
        else:
            raise TypeError("only lists, strings and numbers are supported!")

    def _interactive_pipe_read(self):
        return self._interactive_library_read.readline()

    def calc_static(
        self,
        electronic_steps=100,
        blockSize=8,
        dSpinMoment=1e-8,
        algorithm=None,
        retain_charge_density=False,
        retain_electrostatic_potential=False,
    ):
        """
        Function to setup the hamiltonian to perform static SCF DFT runs

        Args:
            retain_electrostatic_potential:
            retain_charge_density:
            algorithm:
            electronic_steps (int): maximum number of electronic steps, which can be used to achieve convergence
        """
        super(SphinxInteractive, self).calc_static(
            electronic_steps=electronic_steps,
            algorithm=algorithm,
            retain_charge_density=retain_charge_density,
            retain_electrostatic_potential=retain_electrostatic_potential,
        )

    def load_main_group(self):
        main_group = Group()
        if (
            self.server.run_mode.interactive
            or self.server.run_mode.interactive_non_modal
        ):
            commands = Group([
                {
                    "id": '"restart"',
                    "scfDiag":
                        self.get_scf_group(
                            maxSteps=10, keepRhoFixed=True, dEnergy=1.0e-4
                        )
                }, {
                    "id": '"electronicminimization"',
                    "scfDiag": self.get_scf_group(),
                }
            ])
            self.input.sphinx.main.extControl = Group()
            self.input.sphinx.main.extControl.set_group('bornOppenheimer')
            self.input.sphinx.main.extControl.bornOppenheimer = commands
        else:
            super(SphinxInteractive, self).load_main_group()


class SphinxOutput(GenericInteractiveOutput):
    def __init__(self, job):
        super(SphinxOutput, self).__init__(job)

    def check_band_occupancy(self, plot=True):
        """
            Check whether there are still empty bands available.

            args:
                plot (bool): plots occupancy of the last step

            returns:
                True if there are still empty bands
        """
        import matplotlib.pylab as plt
        elec_dict = self._job['output/generic/dft']['n_valence']
        if elec_dict is None:
            raise AssertionError('Number of electrons not parsed')
        n_elec = np.sum([elec_dict[k]
                         for k in self._job.structure.get_chemical_symbols()])
        n_elec = int(np.ceil(n_elec/2))
        bands = self._job['output/generic/dft/bands_occ'][-1]
        bands = bands.reshape(-1, bands.shape[-1])
        max_occ = np.sum(bands>0, axis=-1).max()
        n_bands = bands.shape[-1]
        if plot:
            xticks = np.arange(1, n_bands+1)
            plt.xlabel('Electron number')
            plt.ylabel('Occupancy')
            if n_bands<20:
                plt.xticks(xticks)
            plt.axvline(n_elec, label='#electrons: {}'.format(n_elec))
            plt.axvline(max_occ, color='red',
                label='Max occupancy: {}'.format(max_occ))
            plt.axvline(n_bands, color='green',
                label='Number of bands: {}'.format(n_bands))
            plt.plot(xticks, bands.T, 'x', color='black')
            plt.legend()
        if max_occ < n_bands:
            return True
        else:
            return False


class SphinxInt2(SphinxInteractive):
    def __init__(self, project, job_name):
        warnings.warn("Please use SphinxInt instead of SphinxInt2")
        super(SphinxInt2, self).__init__(project=project, job_name=job_name)
