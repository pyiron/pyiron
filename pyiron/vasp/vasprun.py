# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import os
# import xml.etree.cElementTree as ETree
import numpy as np
from collections import OrderedDict
from pyiron.atomistics.structure.atoms import Atoms
from pyiron.atomistics.structure.periodic_table import PeriodicTable
from pyiron.base.settings.generic import Settings
from pyiron.dft.waves.electronic import ElectronicStructure
import defusedxml.cElementTree as ETree
from defusedxml.ElementTree import ParseError


__author__ = "Sudarsan Surendralal"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - " \
                "Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Sudarsan Surendralal"
__email__ = "surendralal@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


class Vasprun(object):

    """
    This module is used to parse vasprun.xml files and store the data consistent with the pyiron input/output storage
    formats.

    Attributes:

        vasprun_dict (dict): Dictionary containing all information from the calculation parsed from the vasprun.xml
                            file. If you consider a simulation with N atoms and M ionic steps

                'positions' (numpy.ndarray): MxNx3 array containing all the relative positions
                'cell' (numpy.ndarray): Mx3x3 array containing all the size and shape of cells at every iteration point
                'forces' (numpy.ndarray): MxNx3 array containing all the forces in eV/A
                'total_energies' (numpy.ndarray): 1xM array containing all the total energies in eV
    """

    def __init__(self):
        self.vasprun_dict = dict()
        self.root = None

    def from_file(self, filename="vasprun.xml"):
        """
        Parsing vasprun.xml from the working directory

        Args:
            filename (str): Path to the vasprun file
        """
        if not (os.path.isfile(filename)):
            raise AssertionError()
        try:
            self.root = ETree.parse(filename).getroot()
        except ParseError:
            raise VasprunError("The vasprun.xml file is either corrupted or the simulation has failed")

        self.parse_root_to_dict()

    def parse_root_to_dict(self):
        """
        Parses from the main xml root.
        """
        node = self.root
        d = self.vasprun_dict
        d["scf_energies"] = list()
        d["scf_fr_energies"] = list()
        d["scf_0_energies"] = list()
        d["scf_dipole_moments"] = list()
        d["positions"] = list()
        d["cells"] = list()
        d["forces"] = list()
        d["total_energies"] = list()
        d["total_fr_energies"] = list()
        d["total_0_energies"] = list()
        d["stress_tensors"] = list()
        for leaf in node:
            if leaf.tag in ["generator", "incar"]:
                d[leaf.tag] = dict()
                for items in leaf:
                    d[leaf.tag] = self.parse_item_to_dict(items, d[leaf.tag])
            if leaf.tag in ["kpoints"]:
                d[leaf.tag] = dict()
                self.parse_kpoints_to_dict(leaf, d[leaf.tag])
            if leaf.tag in ["atominfo"]:
                d[leaf.tag] = dict()
                self.parse_atom_information_to_dict(leaf, d[leaf.tag])
            if leaf.tag in ["structure"] and leaf.attrib["name"] == "initialpos":
                d["init_structure"] = dict()
                self.parse_structure_to_dict(leaf, d["init_structure"])
            if leaf.tag in ["structure"] and leaf.attrib["name"] == "finalpos":
                d["final_structure"] = dict()
                self.parse_structure_to_dict(leaf, d["final_structure"])
            if leaf.tag in ["calculation"]:
                self.parse_calc_to_dict(leaf, d)
            if leaf.tag in ["parameters"]:
                pass
                self.parse_parameters(leaf, d)
        d["cells"] = np.array(d["cells"])
        d["positions"] = np.array(d["positions"])
        # Check if the parsed coordinates are in absolute/relative coordinates. If absolute, convert to relative
        total_positions = d["positions"].flatten()
        if len(np.argwhere(total_positions > 1)) / len(total_positions) > 0.2:
            pos_new = d["positions"].copy()
            for i, pos in enumerate(pos_new):
                d["positions"][i] = np.dot(pos, np.linalg.inv(d["cells"][i]))
        d["forces"] = np.array(d["forces"])
        d["total_energies"] = np.array(d["total_energies"])
        d["total_fr_energies"] = np.array(d["total_fr_energies"])
        d["total_0_energies"] = np.array(d["total_0_energies"])
        d["scf_energies"] = d["scf_energies"]
        d["scf_dipole_moments"] = d["scf_dipole_moments"]
        d["scf_fr_energies"] = d["scf_fr_energies"]
        d["scf_0_energies"] = d["scf_0_energies"]
        d["stress_tensors"] = d["stress_tensors"]

    def parse_kpoints_to_dict(self, node, d):
        """
        Parses k-points data from a node to a dictionary

        Args:
            node (xml.etree.Element instance): The node to parse
            d (dict): The dictionary to which data is to be parsed
        """
        if not (node.tag == "kpoints"):
            raise AssertionError()
        for leaf in node:
            if leaf.tag == "generation":
                d[leaf.tag] = dict()
                d[leaf.tag]["scheme"] = leaf.attrib["param"]
                if d[leaf.tag]["scheme"] == "listgenerated":
                    line_mode_kpoints = list()
                    for item in leaf:
                        if item.tag == "v":
                            line_mode_kpoints.append(self._parse_vector(item, vec_type=float))
                    d["line_mode_kpoints"] = np.array(line_mode_kpoints)
                else:
                    gen_vec = np.zeros((3, 3))
                    for item in leaf:
                        if item.tag == "v":
                            if item.attrib["name"] in ["divisions"]:
                                d[leaf.tag]["divisions"] = self._parse_vector(item, vec_type=int)
                            if item.attrib["name"] in ["genvec{}".format(i) for i in range(1, 4)]:
                                if item.attrib["name"] == "genvec1":
                                    gen_vec[0, :] = self._parse_vector(item, vec_type=float)
                                if item.attrib["name"] == "genvec2":
                                    gen_vec[1, :] = self._parse_vector(item, vec_type=float)
                                if item.attrib["name"] == "genvec3":
                                    gen_vec[2, :] = self._parse_vector(item, vec_type=float)
                            if item.attrib["name"] in ["shift", "usershift"]:
                                d[leaf.tag][item.attrib["name"]] = self._parse_vector(item, vec_type=float)
                    d[leaf.tag]["genvec"] = np.array(gen_vec)
            if leaf.tag == "varray":
                if leaf.attrib["name"] == "kpointlist":
                    d["kpoint_list"] = self._parse_2d_matrix(leaf, vec_type=float)
                if leaf.attrib["name"] == "weights":
                    d["kpoint_weights"] = self._parse_2d_matrix(leaf, vec_type=float).flatten()

    def parse_atom_information_to_dict(self, node, d):
        """
        Parses atom information from a node to a dictionary

        Args:
            node (xml.etree.Element instance): The node to parse
            d (dict): The dictionary to which data is to be parsed
        """
        if not (node.tag == "atominfo"):
            raise AssertionError()
        species_dict = OrderedDict()
        for leaf in node:
            if leaf.tag == "atoms":
                d["n_atoms"] = self._parse_vector(leaf)[0]
            if leaf.tag == "types":
                d["n_species"] = self._parse_vector(leaf)[0]
            if leaf.tag == "array":
                if leaf.attrib["name"] == "atomtypes":
                    for item in leaf:
                        if item.tag == "set":
                            for sp in item:
                                elements = sp
                                if elements[1].text in species_dict.keys():
                                    pse = PeriodicTable()
                                    count = 1
                                    not_unique = True
                                    species_key = None
                                    while not_unique:
                                        species_key = "_".join([elements[1].text, str(count)])
                                        if species_key not in species_dict.keys():
                                            not_unique = False
                                        else:
                                            count += 1
                                    if species_key is not None:
                                        pse.add_element(elements[1].text, species_key)
                                        special_element = pse.element(species_key)
                                        species_dict[special_element] = dict()
                                        species_dict[special_element]["n_atoms"] = int(elements[0].text)
                                        species_dict[special_element]["valence"] = float(elements[3].text)
                                else:
                                    species_key = elements[1].text
                                    species_dict[species_key] = dict()
                                    species_dict[species_key]["n_atoms"] = int(elements[0].text)
                                    species_dict[species_key]["valence"] = float(elements[3].text)
        d["species_dict"] = species_dict
        species_list = list()
        for key, val in species_dict.items():
            for sp in np.tile([key], species_dict[key]["n_atoms"]):
                species_list.append(clean_character(sp))
        d["species_list"] = species_list

    def parse_fermi_level_to_dict(self, node, d):
        """
        Parses fermi level from a node to a dictionary

        Args:
            node (xml.etree.Element instance): The node to parse
            d (dict): The dictionary to which data is to be parsed
        """
        if not (node.tag == "dos"):
            raise AssertionError()
        for item in node:
            if item.tag == "i":
                self.parse_item_to_dict(item, d)

    def parse_total_dos_to_dict(self, node, d):
        """
        Parses total dos data from a node to a dictionary

        Args:
            node (xml.etree.Element instance): The node to parse
            d (dict): The dictionary to which data is to be parsed
        """
        if not (node.tag == "total"):
            raise AssertionError()
        for item in node:
            if item.tag == "array":
                for ii in item:
                    if ii.tag == "set":
                        spin_dos_energies = list()
                        spin_dos_density = list()
                        spin_dos_idensity = list()
                        for sp in ii:
                            if sp.tag == "set" and "spin" in sp.attrib["comment"]:
                                try:
                                    values = self._parse_2d_matrix(sp, vec_type=float)
                                    dos_energies = values[:, 0]
                                    dos_density = values[:, 1]
                                    dos_idensity = values[:, 2]
                                    spin_dos_energies.append(dos_energies)
                                    spin_dos_density.append(dos_density)
                                    spin_dos_idensity.append(dos_idensity)
                                except ValueError:
                                    pass
                        d["spin_dos_energies"] = np.array(spin_dos_energies)
                        d["spin_dos_density"] = np.array(spin_dos_density)
                        d["spin_dos_idensity"] = np.array(spin_dos_idensity)

    def parse_partial_dos_to_dict(self, node, d):
        """
        Parses partial dos data from a node to a dictionary
    
        Args:
            node (xml.etree.Element instance): The node to parse
            d (dict): The dictionary to which data is to be parsed
        """
        if not (node.tag == "partial"):
            raise AssertionError()
        orbital_dict = dict()
        orbital_index = 0
        for item in node:
            if item.tag == "array":
                for ii in item:
                    if ii.tag == "field":
                        if "energy" not in ii.text:
                            orbital_dict[ii.text.replace(" ", "")] = orbital_index
                        orbital_index += 1
                    if ii.tag == "set":
                        atom_resolved_dos = list()
                        for ion in ii:
                            spin_resolved_dos = list()
                            if ion.tag == "set" and "ion" in ion.attrib["comment"]:
                                for sp in ion:
                                    if sp.tag == "set" and "spin" in sp.attrib["comment"]:
                                        values = self._parse_2d_matrix(sp, vec_type=float)
                                        spin_resolved_dos.append(values[:, 1:])
                            atom_resolved_dos.append(spin_resolved_dos)
                        atom_resolved_dos = np.array(atom_resolved_dos)
                        n_atoms, n_spin, n_densities, n_orbitals = np.shape(atom_resolved_dos)
                        new_grand_dos_matrix = np.zeros((n_spin, n_atoms, n_orbitals, n_densities))
                        for i_spin in range(n_spin):
                            for i_atom in range(n_atoms):
                                for i_orbitals in range(n_orbitals):
                                    new_grand_dos_matrix[i_spin, i_atom, i_orbitals, :] = \
                                        atom_resolved_dos[i_atom, i_spin, :, i_orbitals]
                        d["resolved_dos_matrix"] = new_grand_dos_matrix

    def parse_projected_dos_to_dict(self, node, d):
        """
        Parses partial dos data from a node to a dictionary

        Args:
            node (xml.etree.Element instance): The node to parse
            d (dict): The dictionary to which data is to be parsed
        """
        if not (node.tag == "projected"):
            raise AssertionError()
        orbital_dict = dict()
        orbital_index = 0
        for item in node:
            if item.tag == "array":
                for ii in item:
                    if ii.tag == "field":
                        orbital_dict[ii.text] = orbital_index
                        orbital_index += 1
                    if ii.tag == "set":
                        spin_dos_mat = list()
                        for sp in ii:
                            if sp.tag == "set" and "spin" in sp.attrib["comment"]:
                                kpt_dos_mat = list()
                                for kpt in sp:
                                    band_dos_mat = list()
                                    for band in kpt:
                                        dos_matrix = self._parse_2d_matrix(band, vec_type=float)
                                        band_dos_mat.append(dos_matrix)
                                    kpt_dos_mat.append(band_dos_mat)
                                spin_dos_mat.append(kpt_dos_mat)
                        grand_dos_matrix = np.array(spin_dos_mat)
                        d["grand_dos_matrix"] = grand_dos_matrix
                        d["orbital_dict"] = orbital_dict

    def parse_scf(self, node):
        """
        Parses the total energy and dipole moments for a VASP calculation

        Args:
            node: (xml.etree.Element instance): The node to parse

        Returns:
            d (dict): Dictionary to containing parsed data
        """
        d = dict()
        if not (node.tag == "scstep"):
            raise AssertionError()
        for item in node:
            if item.tag == "energy":
                for i in item:
                    if i.attrib["name"] == "e_wo_entrp":
                        d["scf_energy"] = float(i.text)
                    if i.attrib["name"] == "e_fr_energy":
                        d["scf_fr_energy"] = float(i.text)
                    if i.attrib["name"] == "e_0_energy":
                        d["scf_0_energy"] = float(i.text)
            if item.tag == "dipole":
                for i in item:
                    if i.attrib["name"] == "dipole":
                        d["scf_dipole_moment"] = self._parse_vector(i, vec_type=float)
        return d

    def parse_calc_to_dict(self, node, d):
        """
        Parses ionic step data from a node to a dictionary

        Args:
            node (xml.etree.Element instance): The node to parse
            d (dict): The dictionary to which data is to be parsed
        """
        scf_energies = list()
        scf_fr_energies = list()
        scf_0_energies = list()
        scf_moments = list()
        for item in node:
            if item.tag in ["scstep"]:
                scf_dict = self.parse_scf(item)
                scf_energies.append(scf_dict["scf_energy"])
                scf_fr_energies.append(scf_dict["scf_fr_energy"])
                scf_0_energies.append(scf_dict["scf_0_energy"])
                if "scf_dipole_moment" in scf_dict.keys():
                    scf_moments.append(scf_dict["scf_dipole_moment"])
            if item.tag in ["structure"]:
                struct_dict = dict()
                self.parse_structure_to_dict(item, struct_dict)
                d["positions"].append(struct_dict["positions"])
                d["cells"].append(struct_dict["cell"])
            if item.tag in ["varray"] and item.attrib["name"] == "forces":
                d["forces"].append(self._parse_2d_matrix(item, vec_type=float))
            if item.tag in ["stress"]:
                d["stress_tensors"].append(self._parse_2d_matrix(item, vec_type=float))
            if item.tag == "energy":
                for i in item:
                    if i.attrib["name"] == "e_wo_entrp":
                        d["total_energies"].append(float(i.text))
                    if i.attrib["name"] == "e_fr_energy":
                        d["total_fr_energies"].append(float(i.text))
                    if i.attrib["name"] == "e_0_energy":
                        d["total_0_energies"].append(float(i.text))
                    if i.attrib["name"] == "kinetic":
                        d["kinetic_energies"] = float(i.text)
            if item.tag == "eigenvalues":
                self.parse_eigenvalues_to_dict(item, d)

            if item.tag == "dos":
                self.parse_fermi_level_to_dict(item, d)
                d["efermi"] = float(d["efermi"])
                for i in item:
                    if i.tag == "total":
                        try:
                            self.parse_total_dos_to_dict(i, d)
                        except ValueError:
                            pass
                    if i.tag == "partial":
                        try:
                            self.parse_partial_dos_to_dict(i, d)
                        except ValueError:
                            pass

            if item.tag == "projected":
                self.parse_projected_dos_to_dict(item, d)

        d["scf_energies"].append(scf_energies)
        d["scf_fr_energies"].append(scf_fr_energies)
        d["scf_0_energies"].append(scf_0_energies)
        d["scf_dipole_moments"].append(scf_moments)

    def parse_eigenvalues_to_dict(self, node, d):
        """
        Parses eigenvalue and occupancy data from a node to a dictionary

        Args:
            node (xml.etree.Element instance): The node to parse
            d (dict): The dictionary to which data is to be parsed
        """
        if not (node.tag == "eigenvalues"):
            raise AssertionError()
        grand_eigenvalue_matrix = list()
        grand_occupancy_matrix = list()
        for item in node:
            if item.tag == "array":
                for ii in item:
                    if ii.tag == "set":
                        spin_occ_mat = list()
                        spin_eig_mat = list()
                        for sp in ii:
                            if sp.tag == "set" and "spin" in sp.attrib["comment"]:
                                kpt_eig_mat = list()
                                kpt_occ_mat = list()
                                for kpt in sp:
                                    values = self._parse_2d_matrix(kpt, vec_type=float)
                                    eig_vec = values[:, 0].flatten()
                                    occ_vec = values[:, 1].flatten()
                                    kpt_eig_mat.append(eig_vec)
                                    kpt_occ_mat.append(occ_vec)
                                spin_eig_mat.append(kpt_eig_mat)
                                spin_occ_mat.append(kpt_occ_mat)
                        grand_eigenvalue_matrix = np.array(spin_eig_mat)
                        grand_occupancy_matrix = np.array(spin_occ_mat)
        d["grand_eigenvalue_matrix"] = grand_eigenvalue_matrix
        d["grand_occupancy_matrix"] = grand_occupancy_matrix

    def parse_structure_to_dict(self, node, d):
        """
        Parses structure from a node to a dictionary

        Args:
            node (xml.etree.Element instance): The node to parse
            d (dict): The dictionary to which data is to be parsed
        """
        if not (node.tag == "structure"):
            raise AssertionError()
        for leaf in node:
            if leaf.tag == "crystal":
                for item in leaf:
                    if item.tag == "varray" and item.attrib["name"] == "basis":
                        d["cell"] = self._parse_2d_matrix(item)
            if leaf.tag == "varray" and leaf.attrib["name"] == "positions":
                d["positions"] = self._parse_2d_matrix(leaf)

            if leaf.tag == "varray" and leaf.attrib["name"] == "selective":
                d["selective_dynamics"] = self._parse_2d_matrix(leaf, vec_type=bool)

    @staticmethod
    def parse_item_to_dict(node, d):
        """
        Parses values from an item to a dictionary

        Args:
            node (etree.Element instance): Node to be parsed
            d (dict): The dictionary to which data is to be parsed

        Returns:
            d (dictionary)
        """
        type_dict = {"string": str, "float": float, "int": int, "logical": bool}
        logical_dict = {"T": True, "F": False}
        try:
            if node.attrib["type"] == "logical":
                d[node.attrib["name"]] = logical_dict[node.text.strip()]
            else:
                d[node.attrib["name"]] = type_dict[node.attrib["type"]](node.text)
        except (KeyError, IndexError, ValueError):
            d[node.attrib["name"]] = node.text
        return d

    def parse_parameters(self, node, d):
        """
        Parses parameter data from a node to a dictionary

        Args:
            node (xml.etree.Element instance): The node to parse
            d (dict): The dictionary to which data is to be parsed
        """
        if not (node.tag == "parameters"):
            raise AssertionError()
        self.parse_recursively(node, d, key_name="parameters")

    def parse_recursively(self, node, d, key_name=None):
        """
        Parses recursively from a node to a dictionary

        Args:
            node (xml.etree.Element instance): The node to parse
            d (dict): The dictionary to which data is to be parsed
            key_name (str): Forcefully assign a key name in case it is not present in the xml file
        """
        if not len(node) > 0:
            d[clean_key(node.attrib["name"])] = clean_character(node.text)
            return
        else:
            try:
                if key_name is not None:
                    dict_key = clean_key(key_name)
                else:
                    dict_key = clean_key(node.attrib["name"])
                d[dict_key] = dict()
                for item in node:
                    try:
                        self.parse_item_to_dict(item, d[dict_key][item.attrib["name"]])
                    except (KeyError, ValueError, IndexError):
                        try:
                            self.parse_recursively(item, d[dict_key], item.attrib["name"])
                        except (KeyError, ValueError, IndexError):
                            pass
            except KeyError:
                pass

    def _parse_2d_matrix(self, node, vec_type=float):
        """
        Parses a 2D vector from a node

        Args:
            node (xml.etree.Element instance): The node to parse
            vec_type (type): The type of the vector to be parsed

        Returns:
            numpy.ndarray: The required 2D array/vector
        """
        arr = list()
        for item in node:
            arr.append(self._parse_vector(item, vec_type=vec_type))
        return np.array(arr)

    @staticmethod
    def _parse_vector(node, vec_type=float):
        """
        Parses a 1D vector from a node

        Args:
            node (xml.etree.Element instance): The node to parse
            vec_type (type): The type of the vector to be parsed

        Returns:
            numpy.ndarray: The required 1D array/vector
        """
        txt = node.text
        lst = txt.split()
        logical_dict = {"T": True, "F": False}
        if "type" in node.attrib.keys():
            if node.attrib["type"] == "logical":
                return np.array([logical_dict[l.strip()] for l in lst])
            else:
                return np.array([vec_type(l) for l in lst])
        else:
            return np.array([vec_type(l) for l in lst])

    def get_initial_structure(self):
        """
        Gets the initial structure from the simulation

        Returns:
            pyiron.atomistics.structure.atoms.Atoms: The initial structure

        """
        try:
            el_list = self.vasprun_dict["atominfo"]["species_list"]
            cell = self.vasprun_dict["init_structure"]["cell"]
            positions = self.vasprun_dict["init_structure"]["positions"]
            if len(positions[positions > 1.01]) > 0:
                basis = Atoms(el_list, positions=positions, cell=cell)
            else:
                basis = Atoms(el_list, scaled_positions=positions, cell=cell)
            if "selective_dynamics" in self.vasprun_dict["init_structure"].keys():
                basis.add_tag(selective_dynamics=[True, True, True])
                for i, val in enumerate(self.vasprun_dict["init_structure"]["selective_dynamics"]):
                    basis[i].selective_dynamics = val
            return basis
        except KeyError:
            s = Settings()
            s.logger.warning("The initial structure could not be extracted from vasprun properly")
            return

    def get_final_structure(self):
        """
        Gets the final structure from the simulation

        Returns:
            pyiron.atomistics.structure.atoms.Atoms: The final structure

        """
        try:
            basis = self.get_initial_structure()
            basis.cell = self.vasprun_dict["final_structure"]["cell"]
            positions = self.vasprun_dict["final_structure"]["positions"]
            if len(positions[positions > 1.01]) > 0:
                basis.positions = positions
            else:
                basis.scaled_positions = positions
            return basis
        except (KeyError, AttributeError, ValueError):
            return

    def get_electronic_structure(self):
        """
        Get's the electronic structure from the VASP calculation

        Returns:
            pyiron.atomistics.waves.electronic.ElectronicStructure: The electronic structure object

        """
        es_obj = ElectronicStructure()
        es_obj.kpoint_list = self.vasprun_dict["kpoints"]["kpoint_list"]
        es_obj.kpoint_weights = self.vasprun_dict["kpoints"]["kpoint_weights"]
        es_obj.eigenvalue_matrix = self.vasprun_dict["grand_eigenvalue_matrix"][0, :, :]
        es_obj.occupancy_matrix = self.vasprun_dict["grand_occupancy_matrix"][0, :, :]
        if "grand_dos_matrix" in self.vasprun_dict.keys():
            es_obj.grand_dos_matrix = self.vasprun_dict["grand_dos_matrix"]
        if "efermi" in self.vasprun_dict.keys():
            es_obj.efermi = self.vasprun_dict["efermi"]
        if "spin_dos_energies" in self.vasprun_dict.keys():
            es_obj.dos_energies = self.vasprun_dict["spin_dos_energies"][0]
            es_obj.dos_densities = self.vasprun_dict["spin_dos_density"]
            es_obj.dos_idensities = self.vasprun_dict["spin_dos_idensity"]
        if "resolved_dos_matrix" in self.vasprun_dict.keys():
            es_obj.resolved_densities = self.vasprun_dict["resolved_dos_matrix"]
            es_obj.orbital_dict = self.vasprun_dict["orbital_dict"]
        es_obj.generate_from_matrices()
        return es_obj


def clean_character(a, remove_char=" "):
    """
    Args:
        a (str): String to be cleaned
        remove_char (str): Character to be replaced

    Returns:
        str: The clean string
    """
    if isinstance(a, (str, np.str, np.str_)):
        return a.replace(remove_char, "")
    else:
        return a


def clean_key(a, remove_char=" "):
    """
    Replaces blanck spaces from a string for a dictionary key with "_"
    
    Args:
        a (str): String to be cleaned
        remove_char (str): Character to be replaced

    Returns:
        str: The clean string
    """
    if isinstance(a, (str, np.str, np.str_)):
        return a.replace(remove_char, "_")
    else:
        return a


class VasprunError(ValueError):
    pass
