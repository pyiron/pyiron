# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function

import numpy as np

from pyiron.atomistics.structure.atoms import Atoms
from pyiron.dft.waves.dos import Dos

__author__ = "Sudarsan Surendralal"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH "
    "- Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Sudarsan Surendralal"
__email__ = "surendralal@mpie.de"
__status__ = "development"
__date__ = "Sep 1, 2017"


class ElectronicStructure(object):
    """
    This is a generic module to store electronic structure data in a clean way. Kpoint and Band classes are used to
    store information related to kpoints and bands respectively. Every spin configuration has a set of k-points and
    every k-point has a set of bands associated with it. This is loosely adapted from the `pymatgen electronic_structure
    modules`_. Many of the functions have been substantially modified for pyiron

    .. _pymatgen electronic_structure modules: http://pymatgen.org/pymatgen.electronic_structure.bandstructure.html
    """

    def __init__(self):
        self.kpoints = list()
        self._eigenvalues = list()
        self._occupancies = list()
        self._dos_energies = list()
        self._dos_densities = list()
        self._dos_idensities = list()
        self._eg = None
        self._vbm = None
        self._cbm = None
        self._efermi = None
        self._eigenvalue_matrix = None
        self._occupancy_matrix = None
        self._grand_dos_matrix = None
        self._resolved_densities = None
        self._kpoint_list = list()
        self._kpoint_weights = list()
        self.n_spins = 1
        self._structure = None
        self._orbital_dict = None

    def add_kpoint(self, value, weight):
        """
        Appends a Kpoint() instance to the kpoints attribute

        Args:
            value (list/numpy.ndarray): Value of the k-point in cartesian reciprocal coordinates
            weight (float): The weight of the particular k-point

        """
        kpt_obj = Kpoint()
        kpt_obj.value = value
        kpt_obj.weight = weight
        self.kpoints.append(kpt_obj)

    def get_dos(self, n_bins=100):
        """
        Gives a pyiron.objects.waves.dos.Dos instance

        Args:
            n_bins (int): Number of histogram bins for the dos

        Returns:
            pyiron.objects.waves.dos.Dos: Dos instance

        """
        dos_obj = Dos(n_bins=n_bins, es_obj=self)
        return dos_obj

    @property
    def dos_energies(self):
        """
        numpy.ndarray: A (1xN) vector containing the energies with N grid points
        """
        return self._dos_energies

    @dos_energies.setter
    def dos_energies(self, val):
        self._dos_energies = val

    @property
    def dos_densities(self):
        """
        numpy.ndarray: A (SxN) vector containing the density of states for every spin configuration with S spin
                       configurations and N grid points

        """
        return self._dos_densities

    @dos_densities.setter
    def dos_densities(self, val):
        self._dos_densities = val

    @property
    def dos_idensities(self):
        """
        numpy.ndarray: A (SxN) vector containing the density of states for every spin configuration with S spin
                       configurations and N grid points

        """
        return self._dos_idensities

    @dos_idensities.setter
    def dos_idensities(self, val):
        self._dos_idensities = val

    @property
    def resolved_densities(self):
        """
        numpy.ndarray: A (SxAxOxN) vector containing the density of states for every spin configuration with S spin
                       configurations, A atoms, O orbitals and N grid points. The labels of the orbitals are found on
                       the orbital_dict

        """
        return self._resolved_densities

    @resolved_densities.setter
    def resolved_densities(self, val):
        self._resolved_densities = val

    @property
    def orbital_dict(self):
        """
        dict: A dictionary of the ordering of the orbitals

        Examples:
            >>> self.orbital_dict[0]
            's'

        """
        return self._orbital_dict

    @orbital_dict.setter
    def orbital_dict(self, val):
        self._orbital_dict = val

    @property
    def eigenvalues(self):
        """
        Returns the eigenvalues for each spin value

        numpy.ndarray: Eigenvalues of the bands
        """
        return np.array([val.reshape(-1) for val in self._eigenvalue_matrix])

    @property
    def occupancies(self):
        """
        Returns the occupancies for each spin value

        numpy.ndarray: Occupancies of the bands
        """
        return np.array([val.reshape(-1) for val in self._occupancy_matrix])

    @property
    def eigenvalue_matrix(self):
        """
        numpy.ndarray: A getter function to return the eigenvalue_matrix. The eigenvalue for a given kpoint index i and
                       band index j is given by eigenvalue_matrix[i][j]

        """
        if self._eigenvalue_matrix is None and len(self.kpoints) > 0:
            self._eigenvalue_matrix = np.zeros(
                (len(self.kpoints), len(self.kpoints[0].bands))
            )
            for i, k in enumerate(self.kpoints):
                self._eigenvalue_matrix[i, :] = k.eig_occ_matrix[:, 0]
        return self._eigenvalue_matrix

    @eigenvalue_matrix.setter
    def eigenvalue_matrix(self, val):
        self._eigenvalue_matrix = val

    @property
    def occupancy_matrix(self):
        """
        numpy.ndarray: A getter function to return the occupancy_matrix. The occupancy for a given kpoint index i and
                       band index j is given by occupancy_matrix[i][j]
        """
        if self._occupancy_matrix is None and len(self.kpoints) > 0:
            self._occupancy_matrix = np.zeros(
                (len(self.kpoints), len(self.kpoints[0].bands))
            )
            for i, k in enumerate(self.kpoints):
                self._occupancy_matrix[i, :] = k.eig_occ_matrix[:, 1]
        return self._occupancy_matrix

    @occupancy_matrix.setter
    def occupancy_matrix(self, val):
        self._occupancy_matrix = val

    @property
    def kpoint_list(self):
        """
        list: The list of kpoints in cartesian coordinates
        """
        if len(self._kpoint_list) == 0:
            kpt_lst = list()
            for k in self.kpoints:
                kpt_lst.append(k.value)
            self._kpoint_list = kpt_lst
        return self._kpoint_list

    @kpoint_list.setter
    def kpoint_list(self, val):
        self._kpoint_list = val

    @property
    def kpoint_weights(self):
        """
        list: The weights of the kpoints of the electronic structure in cartesian coordinates
        """
        if len(self._kpoint_weights) == 0:
            kpt_lst = list()
            for k in self.kpoints:
                kpt_lst.append(k.weight)
            self._kpoint_weights = kpt_lst
        return self._kpoint_weights

    @kpoint_weights.setter
    def kpoint_weights(self, val):
        self._kpoint_weights = val

    @property
    def structure(self):
        """
        atomistics.structure.atoms.Atoms: The structure associated with the electronic structure object
        """
        return self._structure

    @structure.setter
    def structure(self, val):
        self._structure = val

    def get_vbm(self, resolution=1e-6):
        """
        Gets the valence band maximum (VBM) of the system for each spin value

        Args:
            resolution (float): An occupancy below this value is considered unoccupied

        Returns:
            dict:
                "value" (float): Absolute energy value of the VBM (eV)
                "kpoint": The Kpoint instance associated with the VBM
                "band": The Band instance associated with the VBM
        """
        vbm_spin_dict = dict()
        n_spins = len(self._eigenvalue_matrix)
        for spin in range(n_spins):
            vbm = None
            vbm_spin_dict[spin] = dict()
            vbm_dict = dict()
            for kpt in self.kpoints:
                for band in kpt.bands[spin]:
                    if band.occupancy > resolution:
                        if vbm is None:
                            vbm = band.eigenvalue
                            vbm_dict["value"] = vbm
                            vbm_dict["kpoint"] = kpt
                            vbm_dict["band"] = band
                        else:
                            if band.eigenvalue > vbm:
                                vbm = band.eigenvalue
                                vbm_dict["value"] = vbm
                                vbm_dict["kpoint"] = kpt
                                vbm_dict["band"] = band
                vbm_spin_dict[spin] = vbm_dict
        return vbm_spin_dict

    def get_cbm(self, resolution=1e-6):
        """
        Gets the conduction band minimum (CBM) of the system for each spin value

        Args:
            resolution (float): An occupancy above this value is considered occupied

        Returns:
            dict:
                "value" (float): Absolute energy value of the CBM (eV)
                 "kpoint": The Kpoint instance associated with the CBM
                 "band": The Band instance associated with the CBM
        """
        cbm_spin_dict = dict()
        n_spins = len(self._eigenvalue_matrix)
        for spin in range(n_spins):
            cbm = None
            cbm_spin_dict[spin] = dict()
            cbm_dict = dict()
            for kpt in self.kpoints:
                for band in kpt.bands[spin]:
                    if band.occupancy <= resolution:
                        if cbm is None:
                            cbm = band.eigenvalue
                            cbm_dict["value"] = cbm
                            cbm_dict["kpoint"] = kpt
                            cbm_dict["band"] = band
                        else:
                            if band.eigenvalue < cbm:
                                cbm = band.eigenvalue
                                cbm_dict["value"] = cbm
                                cbm_dict["kpoint"] = kpt
                                cbm_dict["band"] = band
                cbm_spin_dict[spin] = cbm_dict
        return cbm_spin_dict

    def get_band_gap(self, resolution=1e-6):
        """
        Gets the band gap of the system for each spin value

        Args:
            resolution (float): An occupancy above this value is considered occupied

        Returns:
            dict:
                "band gap" (float): The band gap (eV)
                 "vbm": The dictionary associated with the VBM
                 "cbm": The dictionary associated with the CBM
        """
        gap_dict = {}
        vbm_spin_dict = self.get_vbm(resolution)
        cbm_spin_dict = self.get_cbm(resolution)
        for spin, vbm_dict in vbm_spin_dict.items():
            gap_dict[spin] = dict()
            vbm = vbm_dict["value"]
            cbm = cbm_spin_dict[spin]["value"]
            gap_dict[spin]["band_gap"] = max(0.0, cbm - vbm)
            gap_dict[spin]["vbm"] = vbm_dict
            gap_dict[spin]["cbm"] = cbm_spin_dict[spin]
        return gap_dict

    @property
    def eg(self):
        """
        The band gap for each spin channel

        Returns:
            list: list of band gap values for each spin channel
        """
        self._eg = [val["band_gap"] for val in self.get_band_gap().values()]
        return self._eg

    @eg.setter
    def eg(self, val):
        self._eg = val

    @property
    def vbm(self):
        """
        The Kohn-Sham valence band maximum for each spin channel

        Returns:
            list: list of valence band maximum values for each spin channel
        """
        self._vbm = [val["value"] for val in self.get_vbm().values()]
        return self._vbm

    @vbm.setter
    def vbm(self, val):
        self._vbm = val

    @property
    def cbm(self):
        """
        The Kohn-Sham conduction band minimum for each spin channel

        Returns:
            list: list of conduction band minimum values for each spin channel
        """
        self._cbm = [val["value"] for val in self.get_cbm().values()]
        return self._cbm

    @cbm.setter
    def cbm(self, val):
        self._cbm = val

    @property
    def efermi(self):
        """
        float: The Fermi-level of the system (eV). Please note that in the case of DFT this level is the Kohn-Sham Fermi
               level computed by the DFT code.
        """
        return self._efermi

    @efermi.setter
    def efermi(self, val):
        self._efermi = val

    @property
    def is_metal(self):
        """
        Tells if the given system is metallic or not in each spin channel (up and down respectively).
        The Fermi level crosses bands in the cas of metals but is present in the band gap in the
        case of semi-conductors.

        Returns:

            list: List of boolean values for each spin channel
        """
        if not (self._efermi is not None):
            raise ValueError(
                "e_fermi has to be set before you can determine if the system is metallic or not"
            )
        n_spin, _, n_bands = np.shape(self.eigenvalue_matrix)
        fermi_crossed = [False] * n_spin
        for spin in range(n_spin):
            for i in range(n_bands):
                values = self.eigenvalue_matrix[spin, :, i]
                if (self.efermi < np.max(values)) and (self.efermi >= np.min(values)):
                    fermi_crossed[spin] = True
        return fermi_crossed

    @property
    def grand_dos_matrix(self):
        """
        Getter for the 5 dimensional grand_dos_matrix which gives the contribution of every spin, kpoint, band, atom and
        orbital to the total DOS. For example the dos contribution with spin index s, kpoint k, band b, atom a and
        orbital o is:

        grand_dos_matrix[s, k, b, a, o]

        The grand sum of this matrix would equal 1.0. The spatial, spin, and orbital resolved DOS can be computed using
        this matrix

        Returns:
            numpy.ndarray (5 dimensional)

        """
        if self._grand_dos_matrix is None:
            try:
                n_atoms, n_orbitals = np.shape(
                    self.kpoints[0].bands[0][0].resolved_dos_matrix
                )
            except ValueError:
                return self._grand_dos_matrix
            dimension = (
                self.n_spins,
                len(self.kpoints),
                len(self.kpoints[0].bands),
                n_atoms,
                n_orbitals,
            )
            self._grand_dos_matrix = np.zeros(dimension)
            for spin in range(self.n_spins):
                for i, kpt in enumerate(self.kpoints):
                    for j, band in enumerate(kpt.bands):
                        self._grand_dos_matrix[
                            spin, i, j, :, :
                        ] = band.resolved_dos_matrix
        return self._grand_dos_matrix

    @grand_dos_matrix.setter
    def grand_dos_matrix(self, val):
        """
        Setter for grand_dos_matrix
        """
        self._grand_dos_matrix = val

    def to_hdf(self, hdf, group_name="electronic_structure"):
        """
        Store the object to hdf5 file

        Args:
            hdf: Path to the hdf5 file/group in the file
            group_name: Name of the group under which the attributes are o be stored
        """
        with hdf.open(group_name) as h_es:
            h_es["TYPE"] = str(type(self))
            if self.structure is not None:
                self.structure.to_hdf(h_es)
            h_es["k_points"] = self.kpoint_list
            h_es["k_weights"] = self.kpoint_weights
            h_es["eig_matrix"] = self.eigenvalue_matrix
            h_es["occ_matrix"] = self.occupancy_matrix
            if self.efermi is not None:
                h_es["efermi"] = self.efermi
            with h_es.open("dos") as h_dos:
                h_dos["energies"] = self.dos_energies
                h_dos["tot_densities"] = self.dos_densities
                h_dos["int_densities"] = self.dos_idensities
                if self.grand_dos_matrix is not None:
                    h_dos["grand_dos_matrix"] = self.grand_dos_matrix
                if self.resolved_densities is not None:
                    h_dos["resolved_densities"] = self.resolved_densities

    def from_hdf(self, hdf, group_name="electronic_structure"):
        """
        Retrieve the object from the hdf5 file

        Args:
            hdf: Path to the hdf5 file/group in the file
            group_name: Name of the group under which the attributes are stored
        """
        if "dos" not in hdf[group_name].list_groups():
            self.from_hdf_old(hdf=hdf, group_name=group_name)
        else:
            with hdf.open(group_name) as h_es:
                if "TYPE" not in h_es.list_nodes():
                    h_es["TYPE"] = str(type(self))
                nodes = h_es.list_nodes()
                if self.structure is not None:
                    self.structure.to_hdf(h_es)
                self.kpoint_list = h_es["k_points"]
                self.kpoint_weights = h_es["k_weights"]
                if len(h_es["eig_matrix"].shape) == 2:
                    self.eigenvalue_matrix = np.array([h_es["eig_matrix"]])
                    self.occupancy_matrix = np.array([h_es["occ_matrix"]])
                else:
                    self._eigenvalue_matrix = h_es["eig_matrix"]
                    self._occupancy_matrix = h_es["occ_matrix"]
                self.n_spins = len(self._eigenvalue_matrix)
                if "efermi" in nodes:
                    self.efermi = h_es["efermi"]
                with h_es.open("dos") as h_dos:
                    nodes = h_dos.list_nodes()
                    self.dos_energies = h_dos["energies"]
                    self.dos_densities = h_dos["tot_densities"]
                    self.dos_idensities = h_dos["int_densities"]
                    if "grand_dos_matrix" in nodes:
                        self.grand_dos_matrix = h_dos["grand_dos_matrix"]
                    if "resolved_densities" in nodes:
                        self.resolved_densities = h_dos["resolved_densities"]
            self.generate_from_matrices()

    def to_hdf_old(self, hdf, group_name="electronic_structure"):
        """
        Store the object to hdf5 file

        Args:
            hdf: Path to the hdf5 file/group in the file
            group_name: Name of the group under which the attributes are o be stored
        """
        with hdf.open(group_name) as h_es:
            h_es["TYPE"] = str(type(self))
            if self.structure is not None:
                self.structure.to_hdf(h_es)
            h_es["k_points"] = self.kpoint_list
            h_es["k_point_weights"] = self.kpoint_weights
            h_es["eigenvalue_matrix"] = self.eigenvalue_matrix
            h_es["occupancy_matrix"] = self.occupancy_matrix
            h_es["dos_energies"] = self.dos_energies
            h_es["dos_densities"] = self.dos_densities
            h_es["dos_idensities"] = self.dos_idensities
            if self.efermi is not None:
                h_es["fermi_level"] = self.efermi
            if self.grand_dos_matrix is not None:
                h_es["grand_dos_matrix"] = self.grand_dos_matrix
            if self.resolved_densities is not None:
                h_es["resolved_densities"] = self.resolved_densities

    def from_hdf_old(self, hdf, group_name="electronic_structure"):
        """
        Retrieve the object from the hdf5 file

        Args:
            hdf: Path to the hdf5 file/group in the file
            group_name: Name of the group under which the attributes are stored
        """
        with hdf.open(group_name) as h_es:
            if "structure" in h_es.list_nodes():
                self.structure = Atoms().from_hdf(h_es)
            nodes = h_es.list_nodes()
            self.kpoint_list = h_es["k_points"]
            self.kpoint_weights = h_es["k_point_weights"]
            self.eigenvalue_matrix = np.array([h_es["eigenvalue_matrix"]])
            self.occupancy_matrix = np.array([h_es["occupancy_matrix"]])
            try:
                self.dos_energies = h_es["dos_energies"]
                self.dos_densities = h_es["dos_densities"]
                self.dos_idensities = h_es["dos_idensities"]
            except ValueError:
                pass
            if "fermi_level" in nodes:
                self.efermi = h_es["fermi_level"]
            if "grand_dos_matrix" in nodes:
                self.grand_dos_matrix = h_es["grand_dos_matrix"]
            if "resolved_densities" in nodes:
                self.resolved_densities = h_es["resolved_densities"]
        self.generate_from_matrices()

    def generate_from_matrices(self):
        """
        Generate the Kpoints and Bands from the kpoint lists and sometimes grand_dos_matrix

        """
        for i in range(len(self.kpoint_list)):
            self.add_kpoint(self.kpoint_list[i], self.kpoint_weights[i])
            n_spin, _, length = np.shape(self._eigenvalue_matrix)
            for spin in range(n_spin):
                for j in range(length):
                    val = self.eigenvalue_matrix[spin][i][j]
                    occ = self.occupancy_matrix[spin][i][j]
                    self.kpoints[-1].add_band(eigenvalue=val, occupancy=occ, spin=spin)
                    if self._grand_dos_matrix is not None:
                        self.kpoints[-1].bands[spin][
                            -1
                        ].resolved_dos_matrix = self.grand_dos_matrix[spin, i, j, :, :]

    def get_spin_resolved_dos(self, spin_indices=0):
        """
        Gets the spin resolved DOS

        Args:
            spin_indices (int): The index of the spin for which the DOS is required

        Returns:
            Spin resolved dos (numpy.ndarray instance)

        """
        if not (len(self.dos_energies) > 0):
            raise ValueError("The DOS is not computed/saved for this vasp run")
        return self.dos_densities[spin_indices]

    def get_resolved_dos(self, spin_indices=0, atom_indices=None, orbital_indices=None):
        """
        Get resolved dos based on the specified spin, atom and orbital indices

        Args:
            spin_indices (int/list/numpy.ndarray): spin indices
            atom_indices (int/list/numpy.ndarray): stom indices
            orbital_indices (int/list/numpy.ndarray): orbital indices (based on orbital_dict)

        Returns:
            rdos (numpy.ndarray): Required resolved dos
        """
        if len(self.dos_energies) == 0:
            raise ValueError("The DOS is not computed/saved for this vasp run")
        if self.resolved_densities is None:
            raise ValueError("The resolved DOS is not available for this calculation")
        rdos = None
        if isinstance(spin_indices, (list, np.ndarray)):
            rdos = np.sum(self.resolved_densities[spin_indices], axis=0)
        elif isinstance(spin_indices, int):
            rdos = self.resolved_densities[spin_indices]
        if atom_indices is not None:
            if isinstance(atom_indices, (list, np.ndarray)):
                rdos = np.sum(rdos[atom_indices], axis=0)
            elif isinstance(atom_indices, int):
                rdos = rdos[atom_indices]
        else:
            rdos = np.sum(rdos, axis=0)
        if orbital_indices is not None:
            if isinstance(orbital_indices, (list, np.ndarray)):
                rdos = np.sum(rdos[orbital_indices], axis=0)
            elif isinstance(orbital_indices, int):
                rdos = rdos[orbital_indices]
        else:
            rdos = np.sum(rdos, axis=0)
        return rdos

    def plot_fermi_dirac(self):
        """
        Plots the obtained eigenvalue vs occupation plot

        """
        try:
            import matplotlib.pylab as plt
        except ModuleNotFoundError:
            import matplotlib.pyplot as plt
        for spin, eigenvalues in enumerate(self.eigenvalues):
            arg = np.argsort(eigenvalues)
            plt.plot(eigenvalues[arg], self.occupancies[spin][arg], "-o", label="spin:{}".format(spin), linewidth=2)
        plt.legend()
        plt.axvline(self.efermi, linewidth=2.0, linestyle="dashed", color="black")
        plt.xlabel("Energies (eV)")
        plt.ylabel("Occupancy")
        return plt

    def __del__(self):
        del self.kpoints
        del self._eigenvalues
        del self._occupancies
        del self._eg
        del self._vbm
        del self._cbm
        del self._efermi
        del self._eigenvalue_matrix
        del self._occupancy_matrix
        del self._grand_dos_matrix
        del self._kpoint_list
        del self._kpoint_weights
        del self.n_spins

    def __str__(self):
        output_string = list()
        output_string.append("ElectronicStructure Instance")
        output_string.append("----------------------------")
        output_string.append("Number of spin channels: {}".format(len(self.eigenvalue_matrix)))
        output_string.append("Number of k-points: {}".format(len(self.kpoints)))
        output_string.append("Number of bands: {}".format(len(self.kpoints[0].bands[0])))
        try:
            for spin, is_metal in enumerate(self.is_metal):
                if is_metal:
                    output_string.append("spin {}:".format(spin) + " Is a metal: {}".format(is_metal))
                else:
                    output_string.append("spin {}:".format(spin) + " Is a metal: "
                                                                   "{}".format(is_metal) + " Band gap (ev) "
                                                                                           "{}".format(self.eg[spin]))
        except ValueError:
            pass
        return "\n".join(output_string)

    def __repr__(self):
        return self.__str__()


class Kpoint(object):

    """
    All data related to a single k-point is stored in this module


    Attributes:

        bands (dict): Dict of pyiron.objects.waves.settings.Band objects for each spin channel
        value (float): Value of the k-point
        weight (float): Weight of the k-point used in integration of quantities
        eig_occ_matrix (numpy.ndarray): A Nx2 matrix with the first column with eigenvalues and the second with
                                    occupancies of every band. N being the number of bands assoiated with the k-point
    """

    def __init__(self):
        self._value = None
        self._weight = None
        self.bands = dict()
        self.is_relative = False

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, val):
        self._value = val

    @property
    def weight(self):
        return self._weight

    @weight.setter
    def weight(self, val):
        self._weight = val

    def add_band(self, eigenvalue, occupancy, spin):
        """
        Add a pyiron.objects.waves.core.Band instance

        Args:
            eigenvalue (float): The eigenvalue associated with the Band instance
            occupancy (flaot): The occupancy associated with the Band instance
            spin (int): Spin channel
        """
        band_obj = Band()
        band_obj.eigenvalue = eigenvalue
        band_obj.occupancy = occupancy
        if spin not in self.bands.keys():
            self.bands[spin] = list()
        self.bands[spin].append(band_obj)

    @property
    def eig_occ_matrix(self):
        eig_occ_list = list()
        for bands in self.bands.values():
            eig_occ_list.append([[b.eigenvalue, b.occupancy] for b in bands])
        return np.array(eig_occ_list)


class Band(object):
    """
    All data related to a single band for every k-point is stored in this module
    """

    def __init__(self):
        self._eigenvalue = None
        self._occupancy = None
        self._resolved_dos_matrix = None

    @property
    def eigenvalue(self):
        """
        float: The eigenvalue of a given band at a given k-point
        """
        return self._eigenvalue

    @eigenvalue.setter
    def eigenvalue(self, val):
        self._eigenvalue = val

    @property
    def occupancy(self):
        """
        float: The occupancy of a given band at a given k-point
        """
        return self._occupancy

    @occupancy.setter
    def occupancy(self, val):
        self._occupancy = val

    @property
    def resolved_dos_matrix(self):
        """
        numpy.ndarray instance: 2D matrix with n rows and m columns; n being the unmber of
                                atoms and m being the number of orbitals
        """
        return self._resolved_dos_matrix

    @resolved_dos_matrix.setter
    def resolved_dos_matrix(self, val):
        self._resolved_dos_matrix = val
