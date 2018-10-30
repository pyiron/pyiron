# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
import warnings

__author__ = "Sudarsan Surendralal"
__copyright__ = "Copyright 2017, Max-Planck-Institut für Eisenforschung GmbH " \
                "- Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Sudarsan Surendralal"
__email__ = "surendralal@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"

KBAR_TO_EVA = 6.241509125883258e-4


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
        energies = self.get_total_energies(filename)
        energies_int = self.get_energy_without_entropy(filename)
        energies_zero = self.get_energy_sigma_0(filename)
        scf_energies = self.get_all_total_energies(filename)
        forces = self.get_forces(filename)
        positions = self.get_positions(filename)
        cells = self.get_cells(filename)
        steps = self.get_steps(filename)
        temperatures = self.get_temperatures(filename)
        time = self.get_time(filename)
        fermi_level = self.get_fermi_level(filename)
        scf_moments = self.get_dipole_moments(filename)
        kin_energy_error = self.get_kinetic_energy_error(filename)
        stresses = self.get_stresses(filename, si_unit=False)
        n_elect = self.get_nelect(filename)
        try:
            irreducible_kpoints = self.get_irreducible_kpoints(filename)
        except ValueError:
            print('irreducible kpoints not parsed !')
            irreducible_kpoints = None
        magnetization, final_magmom_lst = self.get_magnetization(filename)
        broyden_mixing = self.get_broyden_mixing_mesh(filename)

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
            hdf (pyiron_base.objects.generic.hdfio.FileHDFio): HDF5 group or file
            group_name (str): Name of the HDF5 group
        """
        with hdf.open(group_name) as hdf5_output:
            for key in self.parse_dict.keys():
                hdf5_output[key] = self.parse_dict[key]

    def to_hdf_minimal(self, hdf, group_name="outcar"):
        """
        Store minimal output in an HDF5 file (output unique to OUTCAR)

        Args:
            hdf (pyiron_base.objects.generic.hdfio.FileHDFio): HDF5 group or file
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
            hdf (pyiron_base.objects.generic.hdfio.FileHDFio): HDF5 group or file
            group_name (str): Name of the HDF5 group
        """
        with hdf.open(group_name) as hdf5_output:
            for key in hdf5_output.list_nodes():
                self.parse_dict[key] = hdf5_output[key]

    @staticmethod
    def get_positions_and_forces(filename="OUTCAR"):
        """
        Gets the forces and positions for every ionic step from the OUTCAR file

        Args:
            filename (str): Filename of the OUTCAR file to parse

        Returns:
            [positions, forces] (sequence)
            numpy.ndarray: A Nx3xM array of positions in $\AA$
            numpy.ndarray: A Nx3xM array of forces in $eV / \AA$

            where N is the number of atoms and M is the number of time steps
        """
        positions = []
        forces = []
        trigger_indices = []
        trigger = "TOTAL-FORCE (eV/Angst)"
        nion_trigger = "NIONS ="
        n_atoms = None
        with open(filename, 'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                line = line.strip()
                if trigger in line:
                    trigger_indices.append(i)
                if n_atoms is None:
                    if nion_trigger in line:
                        n_atoms = int(line.split(nion_trigger)[-1])
        for j in trigger_indices:
            pos = []
            force = []
            for line in lines[j + 2: j + n_atoms + 2]:
                line = line.strip()
                line = _clean_line(line)
                pos.append([float(l) for l in line.split()[0:3]])
                force.append([float(l) for l in line.split()[3:]])
            forces.append(force)
            positions.append(pos)
        return np.array(positions), np.array(forces)

    @staticmethod
    def get_positions(filename="OUTCAR"):

        """
        Gets the positions for every ionic step from the OUTCAR file

        Args:
            filename (str): Filename of the OUTCAR file to parse

        Returns:
            numpy.ndarray: A Nx3xM array of positions in $\AA$

            where N is the number of atoms and M is the number of time steps
        """
        positions = []
        trigger_indices = []
        trigger = "TOTAL-FORCE (eV/Angst)"
        nion_trigger = "NIONS ="
        n_atoms = None
        with open(filename, 'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                line = line.strip()
                if trigger in line:
                    trigger_indices.append(i)
                if n_atoms is None:
                    if nion_trigger in line:
                        n_atoms = int(line.split(nion_trigger)[-1])
        for j in trigger_indices:
            pos = []
            for line in lines[j + 2: j + n_atoms + 2]:
                line = line.strip()
                line = _clean_line(line)
                pos.append([float(l) for l in line.split()[0:3]])
            positions.append(pos)
        return np.array(positions)

    @staticmethod
    def get_forces(filename="OUTCAR"):
        """
        Gets the forces for every ionic step from the OUTCAR file

        Args:
            filename (str): Filename of the OUTCAR file to parse

        Returns:

            numpy.ndarray: A Nx3xM array of forces in $eV / \AA$

            where N is the number of atoms and M is the number of time steps
        """
        forces = []
        trigger_indices = []
        trigger = "TOTAL-FORCE (eV/Angst)"
        nion_trigger = "NIONS ="
        n_atoms = None
        with open(filename, 'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                line = line.strip()
                if trigger in line:
                    trigger_indices.append(i)
                if n_atoms is None:
                    if nion_trigger in line:
                        n_atoms = int(line.split(nion_trigger)[-1])
        for j in trigger_indices:
            force = []
            for line in lines[j + 2: j + n_atoms + 2]:
                line = line.strip()
                line = _clean_line(line)
                force.append([float(l) for l in line.split()[3:]])
            forces.append(force)
        return np.array(forces)

    @staticmethod
    def get_cells(filename="OUTCAR"):
        """
        Gets the cell size and shape for every ionic step from the OUTCAR file

        Args:
            filename (str): Filename of the OUTCAR file to parse

        Returns:
            numpy.ndarray: A 3x3xM array of the cell shape in $\AA$

            where M is the number of time steps
        """
        cells = []
        trigger_indices = []
        trigger = "VOLUME and BASIS-vectors are now :"
        with open(filename, 'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                line = line.strip()
                if trigger in line:
                    trigger_indices.append(i)
        for j in trigger_indices:
            cell = []
            for line in lines[j + 5: j + 8]:
                line = line.strip()
                line = _clean_line(line)
                cell.append([float(l) for l in line.split()[0:3]])
            cells.append(cell)
        return np.array(cells)

    @staticmethod
    def get_stresses(filename="OUTCAR", si_unit=True):
        """

        Args:
            filename (str): Input filename
            si_unit (bool): True SI units are used

        Returns:
            numpy.ndarray: An array of stress values

        """
        trigger = "FORCE on cell =-STRESS in cart. coord.  units (eV):"
        pullay_stress_lst = []
        trigger_indices = []
        with open(filename, 'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                line = line.strip()
                if trigger in line:
                    trigger_indices.append(i)
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
    def get_irreducible_kpoints(filename="OUTCAR", reciprocal=True, weight=True, planewaves=True):
        """
        Function to extract the irreducible kpoints from the OUTCAR file

        Args:
            filename (str): Filename of the OUTCAR file to parse
            reciprocal (bool): Get either the reciprocal or the cartesian coordinates
            weight (bool): Get the weight assigned to the irreducible kpoints
            planewaves (bool): Get the planewaves assigned to the irreducible kpoints

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
        with open(filename, 'r') as f:
            lines = f.readlines()
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
    def get_total_energies(filename="OUTCAR"):
        """
        Gets the total energy for every ionic step from the OUTCAR file

        Args:
            filename (str): Filename of the OUTCAR file to parse

        Returns:
            numpy.ndarray: A 1xM array of the total energies in $eV$

            where M is the number of time steps
        """
        energies = []
        trigger_indices = []
        trigger = "FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)"
        with open(filename, 'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                line = line.strip()
                if trigger in line:
                    trigger_indices.append(i)
        for j in trigger_indices:
            line = lines[j + 2].strip()
            line = _clean_line(line)
            energies.append(float(line.split()[-2]))
        return np.array(energies)

    @staticmethod
    def get_energy_without_entropy(filename="OUTCAR"):
        """
        Gets the total energy for every ionic step from the OUTCAR file

        Args:
            filename (str): Filename of the OUTCAR file to parse

        Returns:
            numpy.ndarray: A 1xM array of the total energies in $eV$

            where M is the number of time steps
        """
        energies = []
        trigger_indices = []
        trigger = "FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)"
        with open(filename, 'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                line = line.strip()
                if trigger in line:
                    trigger_indices.append(i)
        for j in trigger_indices:
            line = lines[j + 4].strip()
            line = _clean_line(line)
            energies.append(float(line.split()[3]))
        return np.array(energies)

    @staticmethod
    def get_energy_sigma_0(filename="OUTCAR"):
        """
        Gets the total energy for every ionic step from the OUTCAR file

        Args:
            filename (str): Filename of the OUTCAR file to parse

        Returns:
            numpy.ndarray: A 1xM array of the total energies in $eV$

            where M is the number of time steps
        """
        energies = []
        trigger_indices = []
        trigger = "FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)"
        with open(filename, 'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                line = line.strip()
                if trigger in line:
                    trigger_indices.append(i)
        for j in trigger_indices:
            line = lines[j + 4].strip()
            line = _clean_line(line)
            energies.append(float(line.split()[-1]))
        return np.array(energies)

    @staticmethod
    def get_all_total_energies(filename="OUTCAR"):
        """
        Gets the energy at every electronic step

        Args:
            filename (str): Filename of the OUTCAR file to parse

        Returns:
            list: A list of energie for every electronic step at every ionic step
        """
        ionic_trigger = "FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)"
        electronic_trigger = "free energy    TOTEN  ="
        scf_energies = list()
        with open(filename, 'r') as f:
            lines = f.readlines()
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
    def get_magnetization(filename="OUTCAR"):
        """
        Gets the magnetization

        Args:
            filename (str): Filename of the OUTCAR file to parse

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
        with open(filename, 'r') as f:
            lines = f.readlines()
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
                if n_atoms is None:
                    if nion_trigger in line:
                        n_atoms = int(line.split(nion_trigger)[-1])
                if local_spin_trigger:
                    for ind_dir, direc in enumerate(['x', 'y', 'z']):
                        if 'magnetization ({})'.format(direc) in line:
                            mag_dict[direc].append([float(lines[i + 4 + atom_index].split()[-1])
                                                    for atom_index in range(n_atoms)])
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
    def get_broyden_mixing_mesh(filename="OUTCAR"):
        """
        Gets the Broyden mixing mesh size

        Args:
            filename (str): Filename of the OUTCAR file to parse

        Returns:
            int: Mesh size
        """
        trigger = "gives a total of "
        with open(filename, 'r') as f:
            for line in f.readlines():
                if trigger in line:
                    return int(line.split()[4])

    @staticmethod
    def get_temperatures(filename="OUTCAR"):
        """
        Gets the temperature at each ionic step (applicable for MD)

        Args:
            filename (str): Filename of the OUTCAR file to parse

        Returns:
            numpy.ndarray: An array of temperatures in Kelvin
        """
        temperatures = []
        trigger_indices = []
        trigger = "kin. lattice  EKIN_LAT= "
        with open(filename, 'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                line = line.strip()
                if trigger in line:
                    trigger_indices.append(i)
        if len(trigger_indices) > 0:
            for j in trigger_indices:
                line = lines[j].strip()
                line = _clean_line(line)
                temperatures.append(float(line.split()[-2]))
        else:
            trigger = "FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)"
            for i, line in enumerate(lines):
                line = line.strip()
                if trigger in line:
                    trigger_indices.append(i)
            temperatures = np.zeros(len(trigger_indices))
        return np.array(temperatures)

    @staticmethod
    def get_steps(filename="OUTCAR"):
        """

        Args:
            filename (str): Filename of the OUTCAR file to parse

        Returns:
            numpy.ndarray: Steps during the simulation
        """
        nblock_trigger = "NBLOCK ="
        trigger = "FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)"
        trigger_indices = list()
        read_nblock = True
        n_block = 1
        with open(filename, 'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                line = line.strip()
                if trigger in line:
                    trigger_indices.append(i)
                if read_nblock is None:
                    if nblock_trigger in line:
                        line = _clean_line(line)
                        n_block = int(line.split(nblock_trigger)[-1])
        return n_block * np.linspace(0, len(trigger_indices))

    def get_time(self, filename="OUTCAR"):
        """
        Time after each simulation step (for MD)

        Args:
            filename (str): Filename of the OUTCAR file to parse

        Returns:
            numpy.ndarray: An array of time values in fs

        """
        potim_trigger = "POTIM  ="
        read_potim = True
        potim = 1.0
        with open(filename, 'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                line = line.strip()
                if read_potim is None:
                    if potim_trigger in line:
                        line = _clean_line(line)
                        potim = float(line.split(potim_trigger)[0])
        return potim * self.get_steps(filename)

    @staticmethod
    def get_kinetic_energy_error(filename="OUTCAR", total=True):
        """
        Get the kinetic energy error

        Args:
            filename (str): Filename of the OUTCAR file to parse
            total (bool): Get either the total correction or the correction per atom

        Returns:
            float: The kinetic energy error in eV
        """
        trigger = "kinetic energy error for atom="
        e_kin_err = None
        n_atoms = None
        nion_trigger = "NIONS ="
        with open(filename, 'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                line = line.strip()
                if trigger in line:
                    e_kin_err = float(line.split()[5])
                if total:
                    if nion_trigger in line:
                        n_atoms = int(line.split(nion_trigger)[-1])
        with open(filename, 'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                line = line.strip()
                if trigger in line:
                    e_kin_err = float(line.split()[5])
        if total and e_kin_err:
            return e_kin_err * n_atoms
        else:
            return e_kin_err

    @staticmethod
    def get_fermi_level(filename="OUTCAR"):
        """
        Getting the Fermi-level (Kohn_Sham) from the OUTCAR file

        Args:
            filename (str): Filename of the OUTCAR file to parse

        Returns:
            float: The Kohn-Sham Fermi level in eV
        """
        trigger = "E-fermi :"
        e_fermi = None
        with open(filename, 'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                line = line.strip()
                if trigger in line:
                    try:
                        e_fermi = float(line.split(trigger)[-1].split()[0])
                    except ValueError:
                        return
        return e_fermi

    @staticmethod
    def get_dipole_moments(filename="OUTCAR"):
        """
        Get the electric dipole moment at every electronic step

        Args:
            filename (str): Filename of the OUTCAR file to parse

        Returns:
            list: A list of dipole moments in (eA) for each electronic step

        """
        moment_trigger = "dipolmoment"
        istep_trigger = "FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)"
        dip_moms = list()
        with open(filename, 'r') as f:
            lines = f.readlines()
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
    def get_nelect(filename="OUTCAR"):
        nelect_trigger = "NELECT"
        with open(filename, 'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                line = line.strip()
                if nelect_trigger in line:
                    return float(line.split()[2])


def _clean_line(line):
    return line.replace("-", " -")
