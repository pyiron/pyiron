# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
import warnings
import scipy.constants
import re

__author__ = "Sudarsan Surendralal"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - " \
                "Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Sudarsan Surendralal"
__email__ = "surendralal@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"

KBAR_TO_EVA = scipy.constants.physical_constants['joule-electron volt relationship'][0] / 1e22


class Outcar(object):
    """
    This module is used to parse VASP OUTCAR files.

    Attributes:

        parse_dict (dict): A dictionary with all the useful quantities parsed from an OUTCAR file after from_file() is
                           executed

    """

    def __init__(self):
        self.parse_dict = dict()

    def from_file(self, filename="OUTCAR"):
        """
        Parse and store relevant quantities from the OUTCAR file into parse_dict.

        Args:
            filename (str): Filename of the OUTCAR file to parse

        """
        with open(filename, 'r') as f:
            lines = f.readlines()
        energies = self.get_total_energies(filename=filename, lines=lines)
        energies_int = self.get_energy_without_entropy(filename=filename, lines=lines)
        energies_zero = self.get_energy_sigma_0(filename=filename, lines=lines)
        scf_energies = self.get_all_total_energies(filename=filename, lines=lines)
        n_atoms = self.get_number_of_atoms(filename=filename, lines=lines)
        forces = self.get_forces(filename=filename, lines=lines, n_atoms=n_atoms)
        positions = self.get_positions(filename=filename, lines=lines, n_atoms=n_atoms)
        cells = self.get_cells(filename=filename, lines=lines)
        steps = self.get_steps(filename=filename, lines=lines)
        temperatures = self.get_temperatures(filename=filename, lines=lines)
        time = self.get_time(filename=filename, lines=lines)
        fermi_level = self.get_fermi_level(filename=filename, lines=lines)
        scf_moments = self.get_dipole_moments(filename=filename, lines=lines)
        kin_energy_error = self.get_kinetic_energy_error(filename=filename, lines=lines)
        stresses = self.get_stresses(filename=filename, si_unit=False, lines=lines)
        n_elect = self.get_nelect(filename=filename, lines=lines)
        try:
            irreducible_kpoints = self.get_irreducible_kpoints(filename=filename, lines=lines)
        except ValueError:
            print('irreducible kpoints not parsed !')
            irreducible_kpoints = None
        magnetization, final_magmom_lst = self.get_magnetization(filename=filename, lines=lines)
        broyden_mixing = self.get_broyden_mixing_mesh(filename=filename, lines=lines)

        self.parse_dict["energies"] = energies
        self.parse_dict["energies_int"] = energies_int
        self.parse_dict["energies_zero"] = energies_zero
        self.parse_dict["scf_energies"] = scf_energies
        self.parse_dict["forces"] = forces
        self.parse_dict["positions"] = positions
        self.parse_dict["cells"] = cells
        self.parse_dict["steps"] = steps
        self.parse_dict["temperatures"] = temperatures
        self.parse_dict["time"] = time
        self.parse_dict["fermi_level"] = fermi_level
        self.parse_dict["scf_dipole_moments"] = scf_moments
        self.parse_dict["kin_energy_error"] = kin_energy_error
        self.parse_dict["stresses"] = stresses
        self.parse_dict["irreducible_kpoints"] = irreducible_kpoints
        self.parse_dict["magnetization"] = magnetization
        self.parse_dict["final_magmoms"] = final_magmom_lst
        self.parse_dict["broyden_mixing"] = broyden_mixing
        self.parse_dict["n_elect"] = n_elect

        try:
            self.parse_dict["pressures"] = np.average(stresses[:, 0:3], axis=1) * KBAR_TO_EVA
        except IndexError:
            self.parse_dict["pressures"] = np.zeros(len(steps))

    def to_hdf(self, hdf, group_name="outcar"):
        """
        Store output in an HDF5 file

        Args:
            hdf (pyiron.base.generic.hdfio.FileHDFio): HDF5 group or file
            group_name (str): Name of the HDF5 group
        """
        with hdf.open(group_name) as hdf5_output:
            for key in self.parse_dict.keys():
                hdf5_output[key] = self.parse_dict[key]

    def to_hdf_minimal(self, hdf, group_name="outcar"):
        """
        Store minimal output in an HDF5 file (output unique to OUTCAR)

        Args:
            hdf (pyiron.base.generic.hdfio.FileHDFio): HDF5 group or file
            group_name (str): Name of the HDF5 group
        """
        unique_quantities = ["kin_energy_error", "broyden_mixing", "stresses", "irreducible_kpoints"]
        with hdf.open(group_name) as hdf5_output:
            for key in self.parse_dict.keys():
                if key in unique_quantities:
                    hdf5_output[key] = self.parse_dict[key]

    def from_hdf(self, hdf, group_name="outcar"):
        """
        Load output from an HDF5 file

        Args:
            hdf (pyiron.base.generic.hdfio.FileHDFio): HDF5 group or file
            group_name (str): Name of the HDF5 group
        """
        with hdf.open(group_name) as hdf5_output:
            for key in hdf5_output.list_nodes():
                self.parse_dict[key] = hdf5_output[key]

    def get_positions_and_forces(self, filename="OUTCAR", lines=None, n_atoms=None):
        """
        Gets the forces and positions for every ionic step from the OUTCAR file

        Args:
            filename (str): Filename of the OUTCAR file to parse
            lines (list/None): lines read from the file
            n_atoms (int/None): number of ions in OUTCAR

        Returns:
            [positions, forces] (sequence)
            numpy.ndarray: A Nx3xM array of positions in $\AA$
            numpy.ndarray: A Nx3xM array of forces in $eV / \AA$

            where N is the number of atoms and M is the number of time steps
        """
        if n_atoms is None:
            n_atoms = self.get_number_of_atoms(filename=filename, lines=lines)
        trigger_indices, lines = _get_trigger(lines=lines, filename=filename, trigger="TOTAL-FORCE (eV/Angst)")
        return self._get_positions_and_forces_parser(lines=lines, trigger_indices=trigger_indices, n_atoms=n_atoms,
                                                     pos_flag=True, force_flag=True)

    def get_positions(self, filename="OUTCAR", lines=None, n_atoms=None):

        """
        Gets the positions for every ionic step from the OUTCAR file

        Args:
            filename (str): Filename of the OUTCAR file to parse
            lines (list/None): lines read from the file
            n_atoms (int/None): number of ions in OUTCAR

        Returns:
            numpy.ndarray: A Nx3xM array of positions in $\AA$

            where N is the number of atoms and M is the number of time steps
        """
        if n_atoms is None:
            n_atoms = self.get_number_of_atoms(filename=filename, lines=lines)
        trigger_indices, lines = _get_trigger(lines=lines, filename=filename, trigger="TOTAL-FORCE (eV/Angst)")
        return self._get_positions_and_forces_parser(lines=lines, trigger_indices=trigger_indices, n_atoms=n_atoms,
                                                     pos_flag=True, force_flag=False)

    def get_forces(self, filename="OUTCAR", lines=None, n_atoms=None):
        """
        Gets the forces for every ionic step from the OUTCAR file

        Args:
            filename (str): Filename of the OUTCAR file to parse
            lines (list/None): lines read from the file
            n_atoms (int/None): number of ions in OUTCAR

        Returns:

            numpy.ndarray: A Nx3xM array of forces in $eV / \AA$

            where N is the number of atoms and M is the number of time steps
        """
        if n_atoms is None:
            n_atoms = self.get_number_of_atoms(filename=filename, lines=lines)
        trigger_indices, lines = _get_trigger(lines=lines, filename=filename, trigger="TOTAL-FORCE (eV/Angst)")
        return self._get_positions_and_forces_parser(lines=lines, trigger_indices=trigger_indices, n_atoms=n_atoms,
                                                     pos_flag=False, force_flag=True)

    def get_cells(self, filename="OUTCAR", lines=None):
        """
        Gets the cell size and shape for every ionic step from the OUTCAR file

        Args:
            filename (str): Filename of the OUTCAR file to parse
            lines (list/None): lines read from the file

        Returns:
            numpy.ndarray: A 3x3xM array of the cell shape in $\AA$

            where M is the number of time steps
        """
        trigger_indices, lines = _get_trigger(lines=lines, filename=filename, trigger="VOLUME and BASIS-vectors are now :")
        return self._get_cells_praser(lines=lines, trigger_indices=trigger_indices)

    @staticmethod
    def get_stresses(filename="OUTCAR", lines=None, si_unit=True):
        """

        Args:
            filename (str): Input filename
            lines (list/None): lines read from the file
            si_unit (bool): True SI units are used

        Returns:
            numpy.ndarray: An array of stress values

        """
        trigger_indices, lines = _get_trigger(lines=lines,
                                              filename=filename,
                                              trigger="FORCE on cell =-STRESS in cart. coord.  units (eV):")
        pullay_stress_lst = []
        for j in trigger_indices:
            try:
                if si_unit:
                    pullay_stress_lst.append([float(l) for l in lines[j + 13].split()[1:7]])
                else:
                    pullay_stress_lst.append([float(l) for l in lines[j + 14].split()[2:8]])
            except ValueError:
                if si_unit:
                    pullay_stress_lst.append([float('NaN')] * 6)
                else:
                    pullay_stress_lst.append([float('NaN')] * 6)
        return np.array(pullay_stress_lst)

    @staticmethod
    def get_irreducible_kpoints(filename="OUTCAR", reciprocal=True, weight=True, planewaves=True, lines=None):
        """
        Function to extract the irreducible kpoints from the OUTCAR file

        Args:
            filename (str): Filename of the OUTCAR file to parse
            reciprocal (bool): Get either the reciprocal or the cartesian coordinates
            weight (bool): Get the weight assigned to the irreducible kpoints
            planewaves (bool): Get the planewaves assigned to the irreducible kpoints
            lines (list/None): lines read from the file

        Returns:
            numpy.ndarray: An array of k-points
        """
        kpoint_lst = []
        weight_lst = []
        planewaves_lst = []
        trigger_number_str = "Subroutine IBZKPT returns following result:"
        trigger_plane_waves_str = "k-point  1 :"
        trigger_number = 0
        trigger_plane_waves = 0
        lines = _get_lines_from_file(filename=filename, lines=lines)
        for i, line in enumerate(lines):
            line = line.strip()
            if trigger_number_str in line:
                trigger_number = int(i)
            elif planewaves:
                if trigger_plane_waves_str in line:
                    trigger_plane_waves = int(i)
        number_irr_kpoints = int(lines[trigger_number + 3].split()[1])
        if reciprocal:
            trigger_start = trigger_number + 7
        else:
            trigger_start = trigger_number + 10 + number_irr_kpoints
        for line in lines[trigger_start: trigger_start + number_irr_kpoints]:
            line = line.strip()
            line = _clean_line(line)
            kpoint_lst.append([float(l) for l in line.split()[0:3]])
            if weight:
                weight_lst.append(float(line.split()[3]))
        if planewaves:
            for line in lines[trigger_plane_waves: trigger_plane_waves + number_irr_kpoints]:
                line = line.strip()
                line = _clean_line(line)
                planewaves_lst.append(float(line.split()[-1]))
        if weight and planewaves:
            return np.array(kpoint_lst), np.array(weight_lst), np.array(planewaves_lst)
        elif weight:
            return np.array(kpoint_lst), np.array(weight_lst)
        elif planewaves:
            return np.array(kpoint_lst), np.array(planewaves_lst)
        else:
            return np.array(kpoint_lst)

    @staticmethod
    def get_total_energies(filename="OUTCAR", lines=None):
        """
        Gets the total energy for every ionic step from the OUTCAR file

        Args:
            filename (str): Filename of the OUTCAR file to parse
            lines (list/None): lines read from the file

        Returns:
            numpy.ndarray: A 1xM array of the total energies in $eV$

            where M is the number of time steps
        """
        def get_total_energies_from_line(line):
            return float(_clean_line(line.strip()).split()[-2])

        trigger_indices, lines = _get_trigger(lines=lines,
                                              filename=filename,
                                              trigger="FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)")
        return np.array([get_total_energies_from_line(lines[j + 2]) for j in trigger_indices])

    @staticmethod
    def get_energy_without_entropy(filename="OUTCAR", lines=None):
        """
        Gets the total energy for every ionic step from the OUTCAR file

        Args:
            filename (str): Filename of the OUTCAR file to parse
            lines (list/None): lines read from the file

        Returns:
            numpy.ndarray: A 1xM array of the total energies in $eV$

            where M is the number of time steps
        """
        def get_energy_without_entropy_from_line(line):
            return float(_clean_line(line.strip()).split()[3])

        trigger_indices, lines = _get_trigger(lines=lines,
                                              filename=filename,
                                              trigger="FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)")
        return np.array([get_energy_without_entropy_from_line(lines[j + 4]) for j in trigger_indices])

    @staticmethod
    def get_energy_sigma_0(filename="OUTCAR", lines=None):
        """
        Gets the total energy for every ionic step from the OUTCAR file

        Args:
            filename (str): Filename of the OUTCAR file to parse
            lines (list/None): lines read from the file

        Returns:
            numpy.ndarray: A 1xM array of the total energies in $eV$

            where M is the number of time steps
        """
        def get_energy_sigma_0_from_line(line):
            return float(_clean_line(line.strip()).split()[-1])

        trigger_indices, lines = _get_trigger(lines=lines,
                                              filename=filename,
                                              trigger="FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)")
        return np.array([get_energy_sigma_0_from_line(lines[j + 4]) for j in trigger_indices])

    @staticmethod
    def get_all_total_energies(filename="OUTCAR", lines=None):
        """
        Gets the energy at every electronic step

        Args:
            filename (str): Filename of the OUTCAR file to parse
            lines (list/None): lines read from the file

        Returns:
            list: A list of energie for every electronic step at every ionic step
        """
        ionic_trigger = "FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)"
        electronic_trigger = "free energy    TOTEN  ="
        scf_energies = list()
        lines = _get_lines_from_file(filename=filename, lines=lines)
        istep_energies = list()
        for i, line in enumerate(lines):
            line = line.strip()
            if ionic_trigger in line:
                scf_energies.append(np.array(istep_energies))
                istep_energies = list()
            if electronic_trigger in line:
                line = _clean_line(line)
                ene = float(line.split()[-2])
                istep_energies.append(ene)
        return scf_energies

    @staticmethod
    def get_magnetization(filename="OUTCAR", lines=None):
        """
        Gets the magnetization

        Args:
            filename (str): Filename of the OUTCAR file to parse
            lines (list/None): lines read from the file

        Returns:
            list: A list with the mgnetization values
        """
        ionic_trigger = "FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)"
        electronic_trigger = "eigenvalue-minimisations"
        nion_trigger = "NIONS ="
        mag_lst = list()
        local_spin_trigger = False
        n_atoms = None
        mag_dict = dict()
        mag_dict['x'] = list()
        mag_dict['y'] = list()
        mag_dict['z'] = list()
        lines = _get_lines_from_file(filename=filename, lines=lines)
        istep_energies = list()
        final_magmom_lst = list()
        for i, line in enumerate(lines):
            line = line.strip()
            if ionic_trigger in line:
                mag_lst.append(np.array(istep_energies))
                istep_energies = list()
            if 'Atomic Wigner-Seitz radii' in line:
                local_spin_trigger = True

            if electronic_trigger in line:
                try:
                    line = lines[i + 2].split('magnetization')[-1]
                    if line != ' \n':
                        spin_str_lst = line.split()
                        spin_str_len = len(spin_str_lst)
                        if spin_str_len == 1:
                            ene = float(line)
                        elif spin_str_len == 3:
                            ene = [float(spin_str_lst[0]), float(spin_str_lst[1]), float(spin_str_lst[2])]
                        else:
                            warnings.warn('Unrecognized spin configuration.')
                            return mag_lst, final_magmom_lst
                        istep_energies.append(ene)
                except ValueError:
                    warnings.warn("Something went wrong in parsing the magnetization")
            if n_atoms is None:
                if nion_trigger in line:
                    n_atoms = int(line.split(nion_trigger)[-1])
            if local_spin_trigger:
                try:
                    for ind_dir, direc in enumerate(['x', 'y', 'z']):
                        if 'magnetization ({})'.format(direc) in line:
                            mag_dict[direc].append([float(lines[i + 4 + atom_index].split()[-1])
                                                    for atom_index in range(n_atoms)])
                except ValueError:
                    warnings.warn("Something went wrong in parsing the magnetic moments")
        if len(mag_dict['x']) > 0:
            if len(mag_dict['y']) == 0:
                final_mag = np.array(mag_dict['x'])
            else:
                n_ionic_steps = np.array(mag_dict['x']).shape[0]
                final_mag = np.abs(np.zeros((n_ionic_steps, n_atoms, 3)))
                final_mag[:, :, 0] = np.array(mag_dict['x'])
                final_mag[:, :, 1] = np.array(mag_dict['y'])
                final_mag[:, :, 2] = np.array(mag_dict['z'])
            final_magmom_lst = final_mag.tolist()
        return mag_lst, final_magmom_lst

    @staticmethod
    def get_broyden_mixing_mesh(filename="OUTCAR", lines=None):
        """
        Gets the Broyden mixing mesh size

        Args:
            filename (str): Filename of the OUTCAR file to parse
            lines (list/None): lines read from the file

        Returns:
            int: Mesh size
        """
        trigger_indices, lines = _get_trigger(lines=lines, filename=filename, trigger="gives a total of ")
        if len(trigger_indices) > 0:
            line_ngx = lines[trigger_indices[0]-2]
        else:
            warnings.warn("Unable to parse the Broyden mixing mesh. Returning 0 instead")
            return 0
        # Exclude all alphabets, and spaces. Then split based on '='
        str_list = re.sub(r'[a-zA-Z]', r'', line_ngx.replace(" ", "").replace("\n", "")).split("=")
        return np.prod([int(val) for val in str_list[1:]])

    @staticmethod
    def get_temperatures(filename="OUTCAR", lines=None):
        """
        Gets the temperature at each ionic step (applicable for MD)

        Args:
            filename (str): Filename of the OUTCAR file to parse
            lines (list/None): lines read from the file

        Returns:
            numpy.ndarray: An array of temperatures in Kelvin
        """
        trigger_indices, lines = _get_trigger(lines=lines, filename=filename, trigger="kin. lattice  EKIN_LAT= ")
        temperatures = []
        if len(trigger_indices) > 0:
            for j in trigger_indices:
                line = lines[j].strip()
                line = _clean_line(line)
                temperatures.append(float(line.split()[-2]))
        else:
            temperatures = np.zeros(len(_get_trigger(lines=lines,
                                                     trigger="FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)",
                                                     return_lines=False)))
        return np.array(temperatures)

    @staticmethod
    def get_steps(filename="OUTCAR", lines=None):
        """

        Args:
            filename (str): Filename of the OUTCAR file to parse
            lines (list/None): lines read from the file

        Returns:
            numpy.ndarray: Steps during the simulation
        """
        nblock_trigger = "NBLOCK ="
        trigger = "FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)"
        trigger_indices = list()
        read_nblock = True
        n_block = 1
        lines = _get_lines_from_file(filename=filename, lines=lines)
        for i, line in enumerate(lines):
            line = line.strip()
            if trigger in line:
                trigger_indices.append(i)
            if read_nblock is None:
                if nblock_trigger in line:
                    line = _clean_line(line)
                    n_block = int(line.split(nblock_trigger)[-1])
        return n_block * np.linspace(0, len(trigger_indices))

    def get_time(self, filename="OUTCAR", lines=None):
        """
        Time after each simulation step (for MD)

        Args:
            filename (str): Filename of the OUTCAR file to parse
            lines (list/None): lines read from the file

        Returns:
            numpy.ndarray: An array of time values in fs

        """
        potim_trigger = "POTIM  ="
        read_potim = True
        potim = 1.0
        lines = _get_lines_from_file(filename=filename, lines=lines)
        for i, line in enumerate(lines):
            line = line.strip()
            if read_potim is None:
                if potim_trigger in line:
                    line = _clean_line(line)
                    potim = float(line.split(potim_trigger)[0])
        return potim * self.get_steps(filename)

    @staticmethod
    def get_kinetic_energy_error(filename="OUTCAR", lines=None):
        """
        Get the kinetic energy error

        Args:
            filename (str): Filename of the OUTCAR file to parse
            lines (list/None): lines read from the file

        Returns:
            float: The kinetic energy error in eV
        """
        trigger = "kinetic energy error for atom="
        e_kin_err = list()
        n_species_list = list()
        nion_trigger = "ions per type ="
        tot_kin_error = 0.0
        lines = _get_lines_from_file(filename=filename, lines=lines)
        for i, line in enumerate(lines):
            line = line.strip()
            if trigger in line:
                e_kin_err.append(float(line.split()[5]))
            if nion_trigger in line:
                n_species_list = [float(val) for val in line.split(nion_trigger)[-1].strip().split()]
        if len(n_species_list) > 0 and len(n_species_list) == len(e_kin_err):
            tot_kin_error = np.sum(np.array(n_species_list) * np.array(e_kin_err))
        return tot_kin_error

    @staticmethod
    def get_fermi_level(filename="OUTCAR", lines=None):
        """
        Getting the Fermi-level (Kohn_Sham) from the OUTCAR file

        Args:
            filename (str): Filename of the OUTCAR file to parse
            lines (list/None): lines read from the file

        Returns:
            float: The Kohn-Sham Fermi level in eV
        """
        trigger = "E-fermi :"
        trigger_indices, lines = _get_trigger(lines=lines, filename=filename, trigger=trigger)
        if len(trigger_indices) != 0:
            try:
                return float(lines[trigger_indices[-1]].split(trigger)[-1].split()[0])
            except ValueError:
                return
        else:
            return

    @staticmethod
    def get_dipole_moments(filename="OUTCAR", lines=None):
        """
        Get the electric dipole moment at every electronic step

        Args:
            filename (str): Filename of the OUTCAR file to parse
            lines (list/None): lines read from the file

        Returns:
            list: A list of dipole moments in (eA) for each electronic step

        """
        moment_trigger = "dipolmoment"
        istep_trigger = "FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)"
        dip_moms = list()
        lines = _get_lines_from_file(filename=filename, lines=lines)
        istep_mom = list()
        for i, line in enumerate(lines):
            line = line.strip()
            if istep_trigger in line:
                dip_moms.append(np.array(istep_mom))
                istep_mom = list()
            if moment_trigger in line:
                line = _clean_line(line)
                mom = np.array([float(val) for val in line.split()[1:4]])
                istep_mom.append(mom)
        return dip_moms

    @staticmethod
    def get_nelect(filename="OUTCAR", lines=None):
        """
        Returns the number of electrons in the simulation

        Args:
            filename (str): OUTCAR filename
            lines (list/None): lines read from the file

        Returns:
            float: The number of electrons in the simulation

        """
        nelect_trigger = "NELECT"
        lines = _get_lines_from_file(filename=filename, lines=lines)
        for i, line in enumerate(lines):
            line = line.strip()
            if nelect_trigger in line:
                return float(line.split()[2])

    @staticmethod
    def get_number_of_atoms(filename="OUTCAR", lines=None):
        """
        Returns the number of ions in the simulation

        Args:
            filename (str): OUTCAR filename
            lines (list/None): lines read from the file

        Returns:
            int: The number of ions in the simulation

        """
        ions_trigger = "NIONS ="
        trigger_indices, lines = _get_trigger(lines=lines, filename=filename, trigger=ions_trigger)
        if len(trigger_indices) != 0:
            return int(lines[trigger_indices[0]].split(ions_trigger)[-1])
        else:
            raise ValueError()

    @staticmethod
    def _get_positions_and_forces_parser(lines, trigger_indices, n_atoms, pos_flag=True, force_flag=True):
        """
        Parser to get the forces and or positions for every ionic step from the OUTCAR file

        Args:
            lines (list): lines read from the file
            trigger_indices (list): list of line indices where the trigger was found.
            n_atoms (int): number of atoms
            pos_flag (bool): parse position
            force_flag (bool): parse forces

        Returns:
            [positions, forces] (sequence)
            numpy.ndarray: A Nx3xM array of positions in $\AA$
            numpy.ndarray: A Nx3xM array of forces in $eV / \AA$

            where N is the number of atoms and M is the number of time steps

        """
        positions = []
        forces = []
        for j in trigger_indices:
            pos = []
            force = []
            for line in lines[j + 2: j + n_atoms + 2]:
                line = line.strip()
                line = _clean_line(line)
                if pos_flag:
                    pos.append([float(l) for l in line.split()[0:3]])
                if force_flag:
                    force.append([float(l) for l in line.split()[3:]])
            forces.append(force)
            positions.append(pos)
        if pos_flag and force_flag:
            return np.array(positions), np.array(forces)
        elif pos_flag:
            return np.array(positions)
        elif force_flag:
            return np.array(forces)

    @staticmethod
    def _get_cells_praser(lines, trigger_indices):
        """
        Parser to get the cell size and shape for every ionic step from the OUTCAR file

        Args:
            lines (list): lines read from the file
            trigger_indices (list): list of line indices where the trigger was found.
            n_atoms (int): number of atoms

        Returns:
            numpy.ndarray: A 3x3xM array of the cell shape in $\AA$

            where M is the number of time steps

        """
        cells = []
        for j in trigger_indices:
            cell = []
            for line in lines[j + 5: j + 8]:
                line = line.strip()
                line = _clean_line(line)
                cell.append([float(l) for l in line.split()[0:3]])
            cells.append(cell)
        return np.array(cells)


def _clean_line(line):
    return line.replace("-", " -")


def _get_trigger(trigger, filename=None, lines=None, return_lines=True):
    """
    Find the lines where a specific trigger appears.

    Args:
        trigger (str): string pattern to search for
        lines (list/None): list of lines
        filename (str/None): file to read lines from

    Returns:
        list: indicies of the lines where the trigger string was found and list of lines
    """
    lines = _get_lines_from_file(filename=filename, lines=lines)
    trigger_indicies = [i for i, line in enumerate(lines) if trigger in line.strip()]
    if return_lines:
        return trigger_indicies, lines
    else:
        return trigger_indicies


def _get_lines_from_file(filename, lines=None):
    """
    If lines is None read the lines from the file with the filename filename.

    Args:
        filename (str): file to read lines from
        lines (list/ None): list of lines

    Returns:
        list: list of lines
    """
    if lines is None:
        with open(filename, 'r') as f:
            lines = f.readlines()
    return lines
