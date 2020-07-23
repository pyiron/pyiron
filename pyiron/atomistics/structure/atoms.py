# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import division, print_function
from ase.atoms import Atoms as ASEAtoms
import ast
from copy import copy
from collections import OrderedDict
from math import cos, sin
import numpy as np
from six import string_types
import warnings
from ase.geometry import get_distances
from matplotlib.colors import rgb2hex
from scipy.interpolate import interp1d
import seekpath
from pyiron.atomistics.structure.atom import Atom
from pyiron.atomistics.structure.sparse_list import SparseArray, SparseList
from pyiron.atomistics.structure.periodic_table import (
    PeriodicTable,
    ChemicalElement,
    ElementColorDictionary,
)
from pyiron.base.settings.generic import Settings
from scipy.spatial import cKDTree, Voronoi
import spglib

__author__ = "Joerg Neugebauer, Sudarsan Surendralal"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Sudarsan Surendralal"
__email__ = "surendralal@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"

s = Settings()


class Atoms(ASEAtoms):
    """
    The Atoms class represents all the information required to describe a structure at the atomic scale. This class is
    written in such a way that is compatible with the `ASE atoms class`_. Some of the functions in this module is based
    on the corresponding implementation in the ASE package

    Args:
        elements (list/numpy.ndarray): List of strings containing the elements or a list of
                            atomistics.structure.periodic_table.ChemicalElement instances
        numbers (list/numpy.ndarray): List of atomic numbers of elements
        symbols (list/numpy.ndarray): List of chemical symbols
        positions (list/numpy.ndarray): List of positions
        scaled_positions (list/numpy.ndarray): List of scaled positions (relative coordinates)
        pbc (list/numpy.ndarray/boolean): Tells if periodic boundary conditions should be applied on the three axes
        cell (list/numpy.ndarray instance): A 3x3 array representing the lattice vectors of the structure

    Note: Only one of elements/symbols or numbers should be assigned during initialization

    Attributes:

        indices (numpy.ndarray): A list of size N which gives the species index of the structure which has N atoms

    .. _ASE atoms class: https://wiki.fysik.dtu.dk/ase/ase/atoms.html

    """

    def __init__(
        self,
        symbols=None,
        positions=None,
        numbers=None,
        tags=None,
        momenta=None,
        masses=None,
        magmoms=None,
        charges=None,
        scaled_positions=None,
        cell=None,
        pbc=None,
        celldisp=None,
        constraint=None,
        calculator=None,
        info=None,
        indices=None,
        elements=None,
        dimension=None,
        species=None,
        **qwargs
    ):
        if symbols is not None:
            if elements is None:
                elements = symbols
            else:
                raise ValueError("Only elements OR symbols should be given.")
        if (
            tags is not None
            or momenta is not None
            or masses is not None
            or charges is not None
            or celldisp is not None
            or constraint is not None
            or calculator is not None
            or info is not None
        ):
            s.logger.debug("Not supported parameter used!")

        self._store_elements = dict()
        self._species_to_index_dict = None
        self.colorLut = ElementColorDictionary().to_lut()
        self._is_scaled = False

        self._species = list()
        self.indices = np.array([])

        self._pse = PeriodicTable()
        self._tag_list = SparseArray()

        el_index_lst = list()
        element_list = None
        if numbers is not None:  # for ASE compatibility
            if not (elements is None):
                raise AssertionError()
            elements = self.numbers_to_elements(numbers)
        if elements is not None:
            el_object_list = None
            if isinstance(elements, str):
                element_list = self.convert_formula(elements)
            elif isinstance(elements, (list, tuple, np.ndarray)):
                if not all([isinstance(el, elements[0].__class__) for el in elements]):
                    object_list = list()
                    for el in elements:
                        if isinstance(el, (str, np.str, np.str_)):
                            object_list.append(self.convert_element(el))
                        if isinstance(el, ChemicalElement):
                            object_list.append(el)
                        if isinstance(el, Atom):
                            object_list.append(el.element)
                        if isinstance(el, (int, np.integer)):
                            # pse = PeriodicTable()
                            object_list.append(self._pse.element(el))
                        el_object_list = object_list

                if len(elements) == 0:
                    element_list = elements
                else:
                    if isinstance(elements[0], (list, tuple, np.ndarray)):
                        elements = np.array(elements).flatten()
                    if isinstance(elements[0], string_types):
                        element_list = elements
                    elif isinstance(elements[0], ChemicalElement):
                        el_object_list = elements
                    elif isinstance(elements[0], Atom):
                        el_object_list = [el.element for el in elements]
                        positions = [el.position for el in elements]
                    elif elements.dtype in [int, np.integer]:
                        el_object_list = self.numbers_to_elements(elements)
                    else:
                        raise ValueError(
                            "Unknown static type for element in list: "
                            + str(type(elements[0]))
                        )

            if el_object_list is None:
                el_object_list = [self.convert_element(el) for el in element_list]

            self.set_species(list(set(el_object_list)))
            # species_to_index_dict = {el: i for i, el in enumerate(self.species)}
            el_index_lst = [self._species_to_index_dict[el] for el in el_object_list]

        elif indices is not None:
            el_index_lst = indices
            self.set_species(species)

        self.indices = np.array(el_index_lst)

        el_lst = [el.Abbreviation if el.Parent is None else el.Parent for el in self.species]
        symbols = np.array([el_lst[el] for el in self.indices])
        super(Atoms, self).__init__(symbols=symbols, positions=positions, numbers=None,
                                    tags=tags, momenta=momenta, masses=masses,
                                    magmoms=magmoms, charges=charges,
                                    scaled_positions=scaled_positions, cell=cell,
                                    pbc=pbc, celldisp=celldisp, constraint=constraint,
                                    calculator=calculator, info=info)
        self._tag_list._length = len(self.positions)
        self.bonds = None
        self.units = {"length": "A", "mass": "u"}
        self._symmetry_dataset = None
        self.set_initial_magnetic_moments(magmoms)
        self._high_symmetry_points = None
        self._high_symmetry_path = None
        self.dimension = dimension
        if len(self.positions) > 0:
            self.dimension = len(self.positions[0])
        else:
            self.dimension = 0

    @property
    def species(self):
        """
        list: A list of atomistics.structure.periodic_table.ChemicalElement instances

        """
        return self._species

    # @species.setter
    def set_species(self, value):
        """
        Setting the species list

        Args:
            value (list): A list atomistics.structure.periodic_table.ChemicalElement instances

        """
        if value is None:
            return
        value = list(value)
        self._species_to_index_dict = {el: i for i, el in enumerate(value)}
        self._species = value[:]
        self._store_elements = {el.Abbreviation: el for el in value}

    @property
    def elements(self):
        """
        numpy.ndarray: A size N list of atomistics.structure.periodic_table.ChemicalElement instances according
                       to the ordering of the atoms in the instance

        """
        return np.array([self.species[el] for el in self.indices])

    def get_high_symmetry_points(self):
        """
        dictionary of high-symmetry points defined for this specific structure.

        Returns:
            dict: high_symmetry_points
        """
        return self._high_symmetry_points

    def _set_high_symmetry_points(self, new_high_symmetry_points):
        """
        Sets new high symmetry points dictionary.

        Args:
            new_high_symmetry_points (dict): new high symmetry points
        """
        if not isinstance(new_high_symmetry_points, dict):
            raise ValueError("has to be dict!")
        self._high_symmetry_points = new_high_symmetry_points

    def add_high_symmetry_points(self, new_points):
        """
        Adds new points to the dict of existing high symmetry points.

        Args:
            new_points (dict): Points to add
        """
        if self.get_high_symmetry_points() is None:
            raise AssertionError("Construct high symmetry points first. Use self.create_line_mode_structure().")
        else:
            self._high_symmetry_points.update(new_points)

    def get_high_symmetry_path(self):
        """
        Path used for band structure calculations

        Returns:
            dict: dict of pathes with start and end points.

        """
        return self._high_symmetry_path

    def _set_high_symmetry_path(self, new_path):
        """
        Sets new list for the high symmetry path used for band structure calculations.

        Args:
            new_path (dict): dictionary of lists of tuples with start and end point.
                E.G. {"my_path": [('Gamma', 'X'), ('X', 'Y')]}
        """
        self._high_symmetry_path = new_path

    def add_high_symmetry_path(self, path):
        """
        Adds a new path to the dictionary of pathes for band structure calculations.

        Args:
            path (dict): dictionary of lists of tuples with start and end point.
                E.G. {"my_path": [('Gamma', 'X'), ('X', 'Y')]}
        """
        if self.get_high_symmetry_path() is None:
            raise AssertionError("Construct high symmetry path first. Use self.create_line_mode_structure().")

        for values_all in path.values():
            for values in values_all:
                if not len(values) == 2:
                    raise ValueError(
                        "'{}' is not a propper trace! It has to contain exactly 2 values! (start and end point)".format(
                            values))
                for v in values:
                    if v not in self.get_high_symmetry_points().keys():
                        raise ValueError("'{}' is not a valid high symmetry point".format(v))

        self._high_symmetry_path.update(path)

    def add_tag(self, *args, **qwargs):
        """
        Add tags to the atoms object.

        Examples:

            For selective dynamics::

            >>> self.add_tag(selective_dynamics=[False, False, False])

        """
        self._tag_list.add_tag(*args, **qwargs)

    # @staticmethod
    def numbers_to_elements(self, numbers):
        """
        Convert atomic numbers in element objects (needed for compatibility with ASE)

        Args:
            numbers (list): List of Element Numbers (as Integers; default in ASE)

        Returns:
            list: A list of elements as needed for pyiron

        """
        # pse = PeriodicTable()  # TODO; extend to internal PSE which can contain additional elements and tags
        atom_number_to_element = {}
        for i_el in set(numbers):
            i_el = int(i_el)
            atom_number_to_element[i_el] = self._pse.element(i_el)
        return [atom_number_to_element[i_el] for i_el in numbers]

    def copy(self):
        """
        Returns a copy of the instance

        Returns:
            pyiron.atomistics.structure.atoms.Atoms: A copy of the instance

        """
        return self.__copy__()

    def to_hdf(self, hdf, group_name="structure"):
        """
        Save the object in a HDF5 file

        Args:
            hdf (pyiron.base.generic.hdfio.FileHDFio): HDF path to which the object is to be saved
            group_name (str):
                Group name with which the object should be stored. This same name should be used to retrieve the object

        """
        # import time
        with hdf.open(group_name) as hdf_structure:
            # time_start = time.time()
            hdf_structure["TYPE"] = str(type(self))
            for el in self.species:
                if isinstance(el.tags, dict):
                    with hdf_structure.open("new_species") as hdf_species:
                        el.to_hdf(hdf_species)
            hdf_structure["species"] = [el.Abbreviation for el in self.species]
            hdf_structure["indices"] = self.indices

            with hdf_structure.open("tags") as hdf_tags:
                for tag in self._tag_list.keys():
                    tag_value = self._tag_list[tag]
                    if isinstance(tag_value, SparseList):
                        tag_value.to_hdf(hdf_tags, tag)
            hdf_structure["units"] = self.units
            hdf_structure["dimension"] = self.dimension

            if self.cell is not None:
                with hdf_structure.open("cell") as hdf_cell:
                    hdf_cell["cell"] = self.cell
                    hdf_cell["pbc"] = self.pbc

            # hdf_structure["coordinates"] = self.positions  # "Atomic coordinates"
            hdf_structure["positions"] = self.positions  # "Atomic coordinates"

            # potentials with explicit bonds (TIP3P, harmonic, etc.)
            if self.bonds is not None:
                hdf_structure["explicit_bonds"] = self.bonds

            # print ('time in atoms.to_hdf: ', time.time() - time_start)

            if self._high_symmetry_points is not None:
                hdf_structure["high_symmetry_points"] = self._high_symmetry_points

            if self._high_symmetry_path is not None:
                hdf_structure["high_symmetry_path"] = self._high_symmetry_path

    def from_hdf(self, hdf, group_name="structure"):
        """
        Retrieve the object from a HDF5 file

        Args:
            hdf (pyiron.base.generic.hdfio.FileHDFio): HDF path to which the object is to be saved
            group_name (str): Group name from which the Atoms object is retreived.

        Returns:
            pyiron_atomistic.structure.atoms.Atoms: The retrieved atoms class

        """
        if "indices" in hdf[group_name].list_nodes():
            with hdf.open(group_name) as hdf_atoms:
                if "new_species" in hdf_atoms.list_groups():
                    with hdf_atoms.open("new_species") as hdf_species:
                        self._pse.from_hdf(hdf_species)

                el_object_list = [
                    self.convert_element(el, self._pse) for el in hdf_atoms["species"]
                ]
                self.indices = hdf_atoms["indices"]
                self._tag_list._length = len(self)

                self.set_species(el_object_list)
                self.bonds = None
                if "explicit_bonds" in hdf_atoms.list_nodes():
                    # print "bonds: "
                    self.bonds = hdf_atoms["explicit_bonds"]

                if "tags" in hdf_atoms.list_groups():
                    with hdf_atoms.open("tags") as hdf_tags:
                        tags = hdf_tags.list_nodes()
                        for tag in tags:
                            # tr_dict = {'0': False, '1': True}
                            if isinstance(hdf_tags[tag], (list, np.ndarray)):
                                my_list = hdf_tags[tag]
                                self._tag_list[tag] = SparseList(
                                    my_list, length=len(self)
                                )

                            else:
                                my_dict = hdf_tags.get_pandas(tag).to_dict()
                                my_dict = {
                                    i: val
                                    for i, val in zip(
                                        my_dict["index"], my_dict["values"]
                                    )
                                }
                                self._tag_list[tag] = SparseList(
                                    my_dict, length=len(self)
                                )

                tr_dict = {1: True, 0: False}
                self.dimension = hdf_atoms["dimension"]
                self.units = hdf_atoms["units"]

                self.cell = None
                if "cell" in hdf_atoms.list_groups():
                    with hdf_atoms.open("cell") as hdf_cell:
                        self.cell = hdf_cell["cell"]
                        self.pbc = hdf_cell["pbc"]

                # Backward compatibility
                position_tag = "positions"
                if position_tag not in hdf_atoms.list_nodes():
                    position_tag = "coordinates"
                if "is_absolute" in hdf_atoms.list_nodes():
                    if not tr_dict[hdf_atoms["is_absolute"]]:
                        self.set_scaled_positions(hdf_atoms[position_tag])
                    else:
                        self.positions = hdf_atoms[position_tag]
                else:
                    self.positions = hdf_atoms[position_tag]

                if "bonds" in hdf_atoms.list_nodes():
                    self.bonds = hdf_atoms["explicit_bonds"]

                self._high_symmetry_points = None
                if "high_symmetry_points" in hdf_atoms.list_nodes():
                    self._high_symmetry_points = hdf_atoms["high_symmetry_points"]

                self._high_symmetry_path = None
                if "high_symmetry_path" in hdf_atoms.list_nodes():
                    self._high_symmetry_path = hdf_atoms["high_symmetry_path"]
                return self

        else:
            return self._from_hdf_old(hdf, group_name)

    def _from_hdf_old(self, hdf, group_name="structure"):
        """
        This function exits merely for the purpose of backward compatibility
        """
        with hdf.open(group_name) as hdf_atoms:
            self._pse = PeriodicTable()
            if "species" in hdf_atoms.list_groups():
                with hdf_atoms.open("species") as hdf_species:
                    self._pse.from_hdf(hdf_species)
            chemical_symbols = np.array(hdf_atoms["elements"], dtype=str)
            el_object_list = [
                self.convert_element(el, self._pse) for el in chemical_symbols
            ]
            self.set_species(list(set(el_object_list)))
            self.indices = [self._species_to_index_dict[el] for el in el_object_list]
            self._tag_list._length = len(self)
            self.bonds = None
            if "explicit_bonds" in hdf_atoms.list_nodes():
                # print "bonds: "
                self.bonds = hdf_atoms["explicit_bonds"]

            if "tags" in hdf_atoms.list_groups():
                with hdf_atoms.open("tags") as hdf_tags:
                    tags = hdf_tags.list_nodes()
                    for tag in tags:
                        # tr_dict = {'0': False, '1': True}
                        if isinstance(hdf_tags[tag], (list, np.ndarray)):
                            my_list = hdf_tags[tag]
                            self._tag_list[tag] = SparseList(my_list, length=len(self))

                        else:
                            my_dict = hdf_tags.get_pandas(tag).to_dict()
                            my_dict = {
                                i: val
                                for i, val in zip(my_dict["index"], my_dict["values"])
                            }
                            self._tag_list[tag] = SparseList(my_dict, length=len(self))

            self.cell = None
            if "cell" in hdf_atoms.list_groups():
                with hdf_atoms.open("cell") as hdf_cell:
                    self.cell = hdf_cell["cell"]
                    self.pbc = hdf_cell["pbc"]

            tr_dict = {1: True, 0: False}
            self.dimension = hdf_atoms["dimension"]
            if "is_absolute" in hdf_atoms and not tr_dict[hdf_atoms["is_absolute"]]:
                self.positions = hdf_atoms["coordinates"]
            else:
                self.set_scaled_positions(hdf_atoms["coordinates"])
            self.units = hdf_atoms["units"]

            if "bonds" in hdf_atoms.list_nodes():
                self.bonds = hdf_atoms["explicit_bonds"]

            self._high_symmetry_points = None
            if "high_symmetry_points" in hdf_atoms.list_nodes():
                self._high_symmetry_points = hdf_atoms["high_symmetry_points"]
            return self

    def center(self, vacuum=None, axis=(0, 1, 2)):
        """
        Center atoms in unit cell.

        Adopted from ASE code (https://wiki.fysik.dtu.dk/ase/_modules/ase/atoms.html#Atoms.center)

        Args:
            vacuum (float): If specified adjust the amount of vacuum when centering. If vacuum=10.0 there will thus be
                            10 Angstrom of vacuum on each side.
            axis (tuple/list): List or turple of integers specifying the axis along which the atoms should be centered

        """

        # Find the orientations of the faces of the unit cell
        c = self.cell
        if c is None:
            c = np.identity(self.dimension)
            self.set_cell(c)

        dirs = np.zeros_like(c)
        for i in range(3):
            dirs[i] = np.cross(c[i - 1], c[i - 2])
            dirs[i] /= np.linalg.norm(dirs[i])  # normalize
            if np.dot(dirs[i], c[i]) < 0.0:
                dirs[i] *= -1

        # Now, decide how much each basis vector should be made longer
        if isinstance(axis, int):
            axes = (axis,)
        else:
            axes = axis
        p = self.positions
        longer = np.zeros(3)
        shift = np.zeros(3)
        for i in axes:
            p0 = np.dot(p, dirs[i]).min()
            p1 = np.dot(p, dirs[i]).max()
            height = np.dot(c[i], dirs[i])
            if vacuum is not None:
                lng = (p1 - p0 + 2 * vacuum) - height
            else:
                lng = 0.0  # Do not change unit cell size!
            top = lng + height - p1
            shf = 0.5 * (top - p0)
            cosphi = np.dot(c[i], dirs[i]) / np.linalg.norm(c[i])
            longer[i] = lng / cosphi
            shift[i] = shf / cosphi

        # Now, do it!
        translation = np.zeros(3)
        for i in axes:
            nowlen = np.sqrt(np.dot(c[i], c[i]))
            cell = self.cell.copy()
            cell[i] *= 1 + longer[i] / nowlen
            self.set_cell(cell)
            translation += shift[i] * c[i] / nowlen
        self.positions += translation
        if self.pbc is None:
            self.pbc = self.dimension * [True]

    def select_index(self, el):
        """
        Returns the indices of a given element in the structure

        Args:
            el (str/atomistics.structures.periodic_table.ChemicalElement/list): Element for which the indices should
                                                                                  be returned
        Returns:
            numpy.ndarray: An array of indices of the atoms of the given element

        """
        if isinstance(el, str):
            return np.where(self.get_chemical_symbols() == el)[0]
        elif isinstance(el, ChemicalElement):
            return np.where([e == el for e in self.get_chemical_elements()])[0]
        if isinstance(el, (list, np.ndarray)):
            if isinstance(el[0], str):
                return np.where(np.isin(self.get_chemical_symbols(), el))[0]
            elif isinstance(el[0], ChemicalElement):
                return np.where([e in el for e in self.get_chemical_elements()])[0]

    def select_parent_index(self, el):
        """
        Returns the indices of a given element in the structure ignoring user defined elements

        Args:
            el (str/atomistics.structures.periodic_table.ChemicalElement): Element for which the indices should
                                                                                  be returned
        Returns:
            numpy.ndarray: An array of indices of the atoms of the given element

        """
        parent_basis = self.get_parent_basis()
        return parent_basis.select_index(el)

    def get_tags(self):
        """
        Returns the keys of the stored tags of the structure

        Returns:
            dict_keys: Keys of the stored tags

        """
        return self._tag_list.keys()

    def convert_element(self, el, pse=None):
        """
        Convert a string or an atom instance into a ChemicalElement instance

        Args:
            el (str/atomistics.structure.atom.Atom): String or atom instance from which the element should
                                                            be generated
            pse (atomistics.structure.periodictable.PeriodicTable): PeriodicTable instance from which the element
                                                                           is generated (optional)

        Returns:

            atomistics.structure.periodictable.ChemicalElement: The required chemical element

        """
        if el in list(self._store_elements.keys()):
            return self._store_elements[el]

        if isinstance(el, string_types):  # as symbol
            element = Atom(el, pse=pse).element
        elif isinstance(el, Atom):
            element = el.element
            el = el.element.Abbreviation
        elif isinstance(el, ChemicalElement):
            element = el
            el = el.Abbreviation
        else:
            raise ValueError("Unknown static type to specify a element")

        self._store_elements[el] = element
        if hasattr(self, "species"):
            if element not in self.species:
                self._species.append(element)
                self.set_species(self._species)
        return element

    def get_chemical_formula(self):
        """
        Returns the chemical formula of structure

        Returns:
            str: The chemical formula as a string

        """
        species = self.get_number_species_atoms()
        formula = ""
        for string_sym, num in species.items():
            if num == 1:
                formula += str(string_sym)
            else:
                formula += str(string_sym) + str(num)
        return formula

    def get_chemical_indices(self):
        """
        Returns the list of chemical indices as ordered in self.species

        Returns:
            numpy.ndarray: A list of chemical indices

        """
        return self.indices

    def get_atomic_numbers(self):
        """
        Returns the atomic numbers of all the atoms in the structure

        Returns:
            numpy.ndarray: A list of atomic numbers

        """
        el_lst = [el.AtomicNumber for el in self.species]
        return np.array([el_lst[el] for el in self.indices])

    def get_chemical_symbols(self):
        """
        Returns the chemical symbols for all the atoms in the structure

        Returns:
            numpy.ndarray: A list of chemical symbols

        """
        el_lst = [el.Abbreviation for el in self.species]
        return np.array([el_lst[el] for el in self.indices])

    def get_parent_symbols(self):
        """
        Returns the chemical symbols for all the atoms in the structure even for user defined elements

        Returns:
            numpy.ndarray: A list of chemical symbols

        """
        sp_parent_list = list()
        for sp in self.species:
            if isinstance(sp.Parent, (float, np.float, type(None))):
                sp_parent_list.append(sp.Abbreviation)
            else:
                sp_parent_list.append(sp.Parent)
        return np.array([sp_parent_list[i] for i in self.indices])

    def get_parent_basis(self):
        """
        Returns the basis with all user defined/special elements as the it's parent

        Returns:
            pyiron.atomistics.structure.atoms.Atoms: Structure without any user defined elements

        """
        parent_basis = copy(self)
        new_species = np.array(parent_basis.species)
        for i, sp in enumerate(new_species):
            if not isinstance(sp.Parent, (float, np.float, type(None))):
                pse = PeriodicTable()
                new_species[i] = pse.element(sp.Parent)
        sym_list = [el.Abbreviation for el in new_species]
        if len(sym_list) != len(np.unique(sym_list)):
            uni, ind, inv_ind = np.unique(
                sym_list, return_index=True, return_inverse=True
            )
            new_species = new_species[ind].copy()
            parent_basis.set_species(list(new_species))
            indices_copy = parent_basis.indices.copy()
            for i, ind_ind in enumerate(inv_ind):
                indices_copy[parent_basis.indices == i] = ind_ind
            parent_basis.indices = indices_copy
            return parent_basis
        parent_basis.set_species(list(new_species))
        return parent_basis

    def get_chemical_elements(self):
        """
        Returns the list of chemical element instances

        Returns:
            numpy.ndarray: A list of chemical element instances

        """
        return self.elements

    def get_number_species_atoms(self):
        """
        Returns a dictionary with the species in the structure and the corresponding count in the structure

        Returns:
            collections.OrderedDict: An ordered dictionary with the species and the corresponding count

        """
        count = OrderedDict()
        # print "sorted: ", sorted(set(self.elements))
        for el in sorted(set(self.get_chemical_symbols())):
            count[el] = 0

        for el in self.get_chemical_symbols():
            count[el] += 1
        return count

    def get_species_symbols(self):
        """
        Returns the symbols of the present species

        Returns:
            numpy.ndarray: List of the symbols of the species

        """
        return np.array(sorted([el.Abbreviation for el in self.species]))

    def get_species_objects(self):
        """


        Returns:

        """
        el_set = self.species
        el_sym_lst = {el.Abbreviation: i for i, el in enumerate(el_set)}
        el_sorted = self.get_species_symbols()
        return [el_set[el_sym_lst[el]] for el in el_sorted]

    def get_number_of_species(self):
        """

        Returns:

        """
        return len(self.species)

    def get_number_of_degrees_of_freedom(self):
        """

        Returns:

        """
        return len(self) * self.dimension

    def get_center_of_mass(self):
        """
        Returns:
            com (float): center of mass in A
        """
        masses = self.get_masses()
        return np.einsum("i,ij->j", masses, self.positions) / np.sum(masses)

    def get_masses(self):
        """

        Returns:

        """
        el_lst = [el.AtomicMass for el in self.species]
        return [el_lst[el] for el in self.indices]

    def get_masses_dof(self):
        """

        Returns:

        """
        dim = self.dimension
        return np.repeat(self.get_masses(), dim)

    def get_volume(self, per_atom=False):
        """

        Args:
            per_atom (bool): True if volume per atom is to be returned

        Returns:
            volume (float): Volume in A**3

        """
        if per_atom:
            return np.abs(np.linalg.det(self.cell)) / len(self)
        else:
            return np.abs(np.linalg.det(self.cell))

    def get_density(self):
        """
        Returns the density in g/cm^3

        Returns:
            float: Density of the structure

        """
        # conv_factor = Ang3_to_cm3/scipi.constants.Avogadro
        # with Ang3_to_cm3 = 1e24
        conv_factor = 1.660539040427164
        return conv_factor * np.sum(self.get_masses()) / self.get_volume()

    def get_number_of_atoms(self):
        """

        Returns:

        """
        # assert(len(self) == np.sum(self.get_number_species_atoms().values()))
        return len(self)

    def set_absolute(self):
        warnings.warn("set_relative is deprecated as of 2020/02/26. It is not guaranteed from v. 0.3", DeprecationWarning)
        if self._is_scaled:
            self._is_scaled = False

    def set_relative(self):
        warnings.warn("set_relative is deprecated as of 2020/02/26. It is not guaranteed from v. 0.3", DeprecationWarning)
        if not self._is_scaled:
            self._is_scaled = True

    def center_coordinates_in_unit_cell(self, origin=0, eps=1e-4):
        """
        Wrap atomic coordinates within the supercell as given by a1, a2., a3

        Args:
            origin (float):  0 to confine between 0 and 1, -0.5 to confine between -0.5 and 0.5
            eps (float): Tolerance to detect atoms at cell edges

        Returns:

            pyiron.atomistics.structure.atoms.Atoms: Wrapped structure

        """
        if any(self.pbc):
            self.set_scaled_positions(
                np.mod(self.get_scaled_positions(wrap=False) + eps, 1) - eps + origin
            )
        return self

    def create_line_mode_structure(self,
                                   with_time_reversal=True,
                                   recipe='hpkot',
                                   threshold=1e-07,
                                   symprec=1e-05,
                                   angle_tolerance=-1.0,
                                   ):
        """
        Uses 'seekpath' to create a new structure with high symmetry points and path for band structure calculations.

        Args:
            with_time_reversal (bool): if False, and the group has no inversion symmetry,
                additional lines are returned as described in the HPKOT paper.
            recipe (str): choose the reference publication that defines the special points and paths.
                Currently, only 'hpkot' is implemented.
            threshold (float): the threshold to use to verify if we are in and edge case
                (e.g., a tetragonal cell, but a==c). For instance, in the tI lattice, if abs(a-c) < threshold,
                a EdgeCaseWarning is issued. Note that depending on the bravais lattice,
                the meaning of the threshold is different (angle, length, …)
            symprec (float): the symmetry precision used internally by SPGLIB
            angle_tolerance (float): the angle_tolerance used internally by SPGLIB

        Returns:
            pyiron.atomistics.structure.atoms.Atoms: new structure
        """
        input_structure = (self.cell, self.get_scaled_positions(), self.indices)
        sp_dict = seekpath.get_path(structure=input_structure,
                                    with_time_reversal=with_time_reversal,
                                    recipe=recipe,
                                    threshold=threshold,
                                    symprec=symprec,
                                    angle_tolerance=angle_tolerance,
                                    )

        original_element_list = [el.Abbreviation for el in self.species]
        element_list = [original_element_list[l] for l in sp_dict["primitive_types"]]
        positions = sp_dict["primitive_positions"]
        pbc = self.pbc
        cell = sp_dict["primitive_lattice"]

        struc_new = Atoms(elements=element_list, scaled_positions=positions, pbc=pbc, cell=cell)

        struc_new._set_high_symmetry_points(sp_dict["point_coords"])
        struc_new._set_high_symmetry_path({"full": sp_dict["path"]})

        return struc_new

    def repeat(self, rep):
        """Create new repeated atoms object.

        The *rep* argument should be a sequence of three positive
        integers like *(2,3,1)* or a single integer (*r*) equivalent
        to *(r,r,r)*."""

        atoms = self.copy()
        atoms *= rep
        return atoms

    def set_repeat(self, vec):
        self *= vec

    def repeat_points(self, points, rep, centered=False):
        """
        Return points with repetition given according to periodic boundary conditions

        Args:
            points (np.ndarray/list): xyz vector or list/array of xyz vectors
            rep (int/list/np.ndarray): Repetition in each direction.
                                       If int is given, the same value is used for
                                       every direction
            centered (bool): Whether the original points should be in the center of
                             repeated points.

        Returns:
            (np.ndarray) repeated points
        """
        n = np.array([rep]).flatten()
        if len(n)==1:
            n = np.tile(n, 3)
        if len(n)!=3:
            raise ValueError('rep must be an integer or a list of 3 integers')
        vector = np.array(points)
        if vector.shape[-1]!=3:
            raise ValueError('points must be an xyz vector or a list/array of xyz vectors')
        if centered and np.mod(n, 2).sum()!=3:
            warnings.warn('When centered, only odd number of repetition should be used')
        v = vector.reshape(-1, 3)
        n_lst = []
        for nn in n:
            if centered:
                n_lst.append(np.arange(nn)-int(nn/2))
            else:
                n_lst.append(np.arange(nn))
        meshgrid = np.meshgrid(n_lst[0], n_lst[1], n_lst[2])
        v_repeated = np.einsum('ni,ij->nj', np.stack(meshgrid, axis=-1).reshape(-1, 3), self.cell)
        v_repeated = v_repeated[:, np.newaxis, :]+v[np.newaxis, :, :]
        return v_repeated.reshape((-1,)+vector.shape)

    def reset_absolute(self, is_absolute):
        raise NotImplementedError("This function was removed!")

    def analyse_ovito_cna_adaptive(self, mode="total"):
        """
        Use Ovito's common neighbor analysis binding.

        Args:
            mode ("total"/"numeric"/"str"): Controls the style and level of detail of the output. (Default is "total", only
                return a summary of the values in the structure.)

        Returns:
            (depends on `mode`)
        """
        from pyiron.atomistics.structure.ovito import analyse_ovito_cna_adaptive
        return analyse_ovito_cna_adaptive(atoms=self, mode=mode)

    def analyse_ovito_centro_symmetry(atoms, num_neighbors=12):
        from pyiron.atomistics.structure.ovito import analyse_ovito_centro_symmetry
        return analyse_ovito_centro_symmetry(atoms, num_neighbors=num_neighbors)

    def analyse_ovito_voronoi_volume(atoms):
        from pyiron.atomistics.structure.ovito import analyse_ovito_voronoi_volume
        return analyse_ovito_voronoi_volume(atoms)

    def analyse_pyscal_steinhardt_parameter(atoms, cutoff=3.5, n_clusters=2, q=[4, 6]):
        from pyiron.atomistics.structure.pyscal import get_steinhardt_parameter_structure
        return get_steinhardt_parameter_structure(structure=atoms, cutoff=cutoff, n_clusters=n_clusters, q=q)

    def analyse_phonopy_equivalent_atoms(atoms):
        from pyiron.atomistics.structure.phonopy import analyse_phonopy_equivalent_atoms

        # warnings.filterwarnings("ignore")
        warnings.warn(
            "analyse_phonopy_equivalent_atoms() is obsolete use get_symmetry()['equivalent_atoms'] instead"
        )
        return analyse_phonopy_equivalent_atoms(atoms)

    @staticmethod
    def _ngl_write_cell(a1, a2, a3, f1=90, f2=90, f3=90):
        """
        Writes a PDB-formatted line to represent the simulation cell.

        Args:
            a1, a2, a3 (float): Lengths of the cell vectors.
            f1, f2, f3 (float): Angles between the cell vectors (which angles exactly?) (in degrees).

        Returns:
            (str): The line defining the cell in PDB format.
        """
        return "CRYST1 {:8.3f} {:8.3f} {:8.3f} {:6.2f} {:6.2f} {:6.2f} P 1\n".format(
            a1, a2, a3, f1, f2, f3
        )

    @staticmethod
    def _ngl_write_atom(
        num,
        species,
        x,
        y,
        z,
        group=None,
        num2=None,
        occupancy=1.0,
        temperature_factor=0.0,
    ):
        """
        Writes a PDB-formatted line to represent an atom.

        Args:
            num (int): Atomic index.
            species (str): Elemental species.
            x, y, z (float): Cartesian coordinates of the atom.
            group (str): A...group name? (Default is None, repeat elemental species.)
            num2 (int): An "alternate" index. (Don't ask me...) (Default is None, repeat first number.)
            occupancy (float): PDB occupancy parameter. (Default is 1.)
            temperature_factor (float): PDB temperature factor parameter. (Default is 0.

        Returns:
            (str): The line defining an atom in PDB format

        Warnings:
            * The [PDB docs](https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html) indicate that
                the xyz coordinates might need to be in some sort of orthogonal basis. If you have weird behaviour,
                this might be a good place to investigate.
        """
        if group is None:
            group = species
        if num2 is None:
            num2 = num
        return "ATOM {:>6} {:>4} {:>4} {:>5} {:10.3f} {:7.3f} {:7.3f} {:5.2f} {:5.2f} {:>11} \n".format(
            num, species, group, num2, x, y, z, occupancy, temperature_factor, species
        )

    def _ngl_write_structure(self, elements, positions, cell):
        """
        Turns structure information into a NGLView-readable protein-database-formatted string.

        Args:
            elements (numpy.ndarray/list): Element symbol for each atom.
            positions (numpy.ndarray/list): Vector of Cartesian atom positions.
            cell (numpy.ndarray/list): Simulation cell Bravais matrix.

        Returns:
            (str): The PDB-formatted representation of the structure.
        """
        from ase.geometry import cell_to_cellpar, cellpar_to_cell
        if cell is None or any(np.max(cell, axis=0) < 1e-2):
            # Define a dummy cell if it doesn't exist (eg. for clusters)
            max_pos = np.max(positions, axis=0)
            max_pos[np.abs(max_pos) < 1e-2] = 10
            cell = np.eye(3) * max_pos
        cellpar = cell_to_cellpar(cell)
        exportedcell = cellpar_to_cell(cellpar)
        rotation = np.linalg.solve(cell, exportedcell)

        pdb_str = self._ngl_write_cell(*cellpar)
        pdb_str += "MODEL     1\n"

        if rotation is not None:
            positions = np.array(positions).dot(rotation)

        for i, p in enumerate(positions):
            pdb_str += self._ngl_write_atom(i, elements[i], *p)

        pdb_str += "ENDMDL \n"
        return pdb_str

    def _atomic_number_to_radius(self, atomic_number, shift=0.2, slope=0.1, scale=1.0):
        """
        Give the atomic radius for plotting, which scales like the root of the atomic number.

        Args:
            atomic_number (int/float): The atomic number.
            shift (float): A constant addition to the radius. (Default is 0.2.)
            slope (float): A multiplier for the root of the atomic number. (Default is 0.1)
            scale (float): How much to rescale the whole thing by.

        Returns:
            (float): The radius. (Not physical, just for visualization!)
        """
        return (shift + slope * np.sqrt(atomic_number)) * scale

    def _add_colorscheme_spacefill(
        self, view, elements, atomic_numbers, particle_size, scheme="element"
    ):
        """
        Set NGLView spacefill parameters according to a color-scheme.

        Args:
            view (NGLWidget): The widget to work on.
            elements (numpy.ndarray/list): Elemental symbols.
            atomic_numbers (numpy.ndarray/list): Integer atomic numbers for determining atomic size.
            particle_size (float): A scale factor for the atomic size.
            scheme (str): The scheme to use. (Default is "element".)

            Possible NGLView color schemes:
              " ", "picking", "random", "uniform", "atomindex", "residueindex",
              "chainindex", "modelindex", "sstruc", "element", "resname", "bfactor",
              "hydrophobicity", "value", "volume", "occupancy"

        Returns:
            (nglview.NGLWidget): The modified widget.
        """
        for elem, num in set(list(zip(elements, atomic_numbers))):
            view.add_spacefill(
                selection="#" + elem,
                radius_type="vdw",
                radius=self._atomic_number_to_radius(num, scale=particle_size),
                color_scheme=scheme,
            )
        return view

    def _add_custom_color_spacefill(self, view, atomic_numbers, particle_size, colors):
        """
        Set NGLView spacefill parameters according to per-atom colors.

        Args:
            view (NGLWidget): The widget to work on.
            atomic_numbers (numpy.ndarray/list): Integer atomic numbers for determining atomic size.
            particle_size (float): A scale factor for the atomic size.
            colors (numpy.ndarray/list): A per-atom list of HTML or hex color codes.

        Returns:
            (nglview.NGLWidget): The modified widget.
        """
        for n, num in enumerate(atomic_numbers):
            view.add_spacefill(
                selection=[n],
                radius_type="vdw",
                radius=self._atomic_number_to_radius(num, scale=particle_size),
                color=colors[n],
            )
        return view

    @staticmethod
    def _scalars_to_hex_colors(scalar_field, start=None, end=None, cmap=None):
        """
        Convert scalar values to hex codes using a colormap.

        Args:
            scalar_field (numpy.ndarray/list): Scalars to convert.
            start (float): Scalar value to map to the bottom of the colormap (values below are clipped). (Default is
                None, use the minimal scalar value.)
            end (float): Scalar value to map to the top of the colormap (values above are clipped).  (Default is
                None, use the maximal scalar value.)
            cmap (matplotlib.cm): The colormap to use. (Default is None, which gives a blue-red divergent map.)

        Returns:
            (list): The corresponding hex codes for each scalar value passed in.
        """
        if start is None:
            start = np.amin(scalar_field)
        if end is None:
            end = np.amax(scalar_field)
        interp = interp1d([start, end], [0, 1])
        remapped_field = interp(
            np.clip(scalar_field, start, end)
        )  # Map field onto [0,1]

        if cmap is None:
            try:
                from seaborn import diverging_palette
            except ImportError:
                print(
                    "The package seaborn needs to be installed for the plot3d() function!"
                )
            cmap = diverging_palette(245, 15, as_cmap=True)  # A nice blue-red palette

        return [
            rgb2hex(cmap(scalar)[:3]) for scalar in remapped_field
        ]  # The slice gets RGB but leaves alpha

    def plot3d(
        self,
        show_cell=True,
        show_axes=True,
        camera="orthographic",
        spacefill=True,
        particle_size=1.0,
        select_atoms=None,
        background="white",
        color_scheme=None,
        colors=None,
        scalar_field=None,
        scalar_start=None,
        scalar_end=None,
        scalar_cmap=None,
        vector_field=None,
        vector_color=None,
        magnetic_moments=False,
        custom_array=None,
        custom_3darray=None,
    ):
        """
        Plot3d relies on NGLView to visualize atomic structures. Here, we construct a string in the "protein database"
        ("pdb") format, then turn it into an NGLView "structure". PDB is a white-space sensitive format, so the
        string snippets are carefully formatted.

        The final widget is returned. If it is assigned to a variable, the visualization is suppressed until that
        variable is evaluated, and in the meantime more NGL operations can be applied to it to modify the visualization.

        Args:
            show_cell (bool): Whether or not to show the frame. (Default is True.)
            show_axes (bool): Whether or not to show xyz axes. (Default is True.)
            camera (str): 'perspective' or 'orthographic'. (Default is 'perspective'.)
            spacefill (bool): Whether to use a space-filling or ball-and-stick representation. (Default is True, use
                space-filling atoms.)
            particle_size (float): Size of the particles. (Default is 1.)
            select_atoms (numpy.ndarray): Indices of atoms to show, either as integers or a boolean array mask.
                (Default is None, show all atoms.)
            background (str): Background color. (Default is 'white'.)
            color_scheme (str): NGLView color scheme to use. (Default is None, color by element.)
            colors (numpy.ndarray): A per-atom array of HTML color names or hex color codes to use for atomic colors.
                (Default is None, use coloring scheme.)
            scalar_field (numpy.ndarray): Color each atom according to the array value (Default is None, use coloring
                scheme.)
            scalar_start (float): The scalar value to be mapped onto the low end of the color map (lower values are
                clipped). (Default is None, use the minimum value in `scalar_field`.)
            scalar_end (float): The scalar value to be mapped onto the high end of the color map (higher values are
                clipped). (Default is None, use the maximum value in `scalar_field`.)
            scalar_cmap (matplotlib.cm): The colormap to use. (Default is None, giving a blue-red divergent map.)
            vector_field (numpy.ndarray): Add vectors (3 values) originating at each atom. (Default is None, no
                vectors.)
            vector_color (numpy.ndarray): Colors for the vectors (only available with vector_field). (Default is None,
                vectors are colored by their direction.)
            magnetic_moments (bool): Plot magnetic moments as 'scalar_field' or 'vector_field'.

            Possible NGLView color schemes:
              " ", "picking", "random", "uniform", "atomindex", "residueindex",
              "chainindex", "modelindex", "sstruc", "element", "resname", "bfactor",
              "hydrophobicity", "value", "volume", "occupancy"

        Returns:
            (nglview.NGLWidget): The NGLView widget itself, which can be operated on further or viewed as-is.

        Warnings:
            * Many features only work with space-filling atoms (e.g. coloring by a scalar field).
            * The colour interpretation of some hex codes is weird, e.g. 'green'.
        """
        try:  # If the graphical packages are not available, the GUI will not work.
            import nglview
        except ImportError:
            raise ImportError(
                "The package nglview needs to be installed for the plot3d() function!"
            )

        if custom_array is not None:
            warnings.warn(
                "custom_array is deprecated. Use scalar_field instead",
                DeprecationWarning,
            )
            scalar_field = custom_array

        if custom_3darray is not None:
            warnings.warn(
                "custom_3darray is deprecated. Use vector_field instead",
                DeprecationWarning,
            )
            vector_field = custom_3darray

        if magnetic_moments is True and hasattr(self, 'spin'):
            if len(self.get_initial_magnetic_moments().shape) == 1:
                scalar_field = self.get_initial_magnetic_moments()
            else:
                vector_field = self.get_initial_magnetic_moments()

        parent_basis = self.get_parent_basis()
        elements = parent_basis.get_chemical_symbols()
        atomic_numbers = parent_basis.get_atomic_numbers()
        positions = self.positions

        # If `select_atoms` was given, visualize only a subset of the `parent_basis`
        if select_atoms is not None:
            select_atoms = np.array(select_atoms, dtype=int)
            elements = elements[select_atoms]
            atomic_numbers = atomic_numbers[select_atoms]
            positions = positions[select_atoms]
            if colors is not None:
                colors = np.array(colors)
                colors = colors[select_atoms]
            if scalar_field is not None:
                scalar_field = np.array(scalar_field)
                scalar_field = scalar_field[select_atoms]
            if vector_field is not None:
                vector_field = np.array(vector_field)
                vector_field = vector_field[select_atoms]
            if vector_color is not None:
                vector_color = np.array(vector_color)
                vector_color = vector_color[select_atoms]

        # Write the nglview protein-database-formatted string
        struct = nglview.TextStructure(
            self._ngl_write_structure(elements, positions, self.cell)
        )

        # Parse the string into the displayable widget
        view = nglview.NGLWidget(struct)

        if spacefill:
            # Color by scheme
            if color_scheme is not None:
                if colors is not None:
                    warnings.warn("`color_scheme` is overriding `colors`")
                if scalar_field is not None:
                    warnings.warn("`color_scheme` is overriding `scalar_field`")
                view = self._add_colorscheme_spacefill(
                    view, elements, atomic_numbers, particle_size, color_scheme
                )
            # Color by per-atom colors
            elif colors is not None:
                if scalar_field is not None:
                    warnings.warn("`colors` is overriding `scalar_field`")
                view = self._add_custom_color_spacefill(
                    view, atomic_numbers, particle_size, colors
                )
            # Color by per-atom scalars
            elif scalar_field is not None:  # Color by per-atom scalars
                colors = self._scalars_to_hex_colors(
                    scalar_field, scalar_start, scalar_end, scalar_cmap
                )
                view = self._add_custom_color_spacefill(
                    view, atomic_numbers, particle_size, colors
                )
            # Color by element
            else:
                view = self._add_colorscheme_spacefill(
                    view, elements, atomic_numbers, particle_size
                )
            view.remove_ball_and_stick()
        else:
            view.add_ball_and_stick()

        if show_cell:
            if parent_basis.cell is not None:
                if all(np.max(parent_basis.cell, axis=0) > 1e-2):
                    view.add_unitcell()

        if vector_color is None and vector_field is not None:
            vector_color = (
                0.5
                * vector_field
                / np.linalg.norm(vector_field, axis=-1)[:, np.newaxis]
                + 0.5
            )
        elif (
            vector_field is not None and vector_field is not None
        ):  # WARNING: There must be a bug here...
            try:
                if vector_color.shape != np.ones((len(self), 3)).shape:
                    vector_color = np.outer(
                        np.ones(len(self)), vector_color / np.linalg.norm(vector_color)
                    )
            except AttributeError:
                vector_color = np.ones((len(self), 3)) * vector_color

        if vector_field is not None:
            for arr, pos, col in zip(vector_field, positions, vector_color):
                view.shape.add_arrow(list(pos), list(pos + arr), list(col), 0.2)

        if show_axes:  # Add axes
            axes_origin = -np.ones(3)
            arrow_radius = 0.1
            text_size = 1
            text_color = [0, 0, 0]
            arrow_names = ["x", "y", "z"]

            for n in [0, 1, 2]:
                start = list(axes_origin)
                shift = np.zeros(3)
                shift[n] = 1
                end = list(start + shift)
                color = list(shift)
                # We cast as list to avoid JSON warnings
                view.shape.add_arrow(start, end, color, arrow_radius)
                view.shape.add_text(end, text_color, text_size, arrow_names[n])

        if camera != "perspective" and camera != "orthographic":
            warnings.warn(
                "Only perspective or orthographic is (likely to be) permitted for camera"
            )

        view.camera = camera
        view.background = background
        return view

    def plot3d_ase(
        self,
        spacefill=True,
        show_cell=True,
        camera="perspective",
        particle_size=0.5,
        background="white",
        color_scheme="element",
        show_axes=True,
    ):
        """
        Possible color schemes:
          " ", "picking", "random", "uniform", "atomindex", "residueindex",
          "chainindex", "modelindex", "sstruc", "element", "resname", "bfactor",
          "hydrophobicity", "value", "volume", "occupancy"
        Returns:
        """
        try:  # If the graphical packages are not available, the GUI will not work.
            import nglview
        except ImportError:
            raise ImportError(
                "The package nglview needs to be installed for the plot3d() function!"
            )
        # Always visualize the parent basis
        parent_basis = self.get_parent_basis()
        view = nglview.show_ase(parent_basis)
        if spacefill:
            view.add_spacefill(
                radius_type="vdw", color_scheme=color_scheme, radius=particle_size
            )
            # view.add_spacefill(radius=1.0)
            view.remove_ball_and_stick()
        else:
            view.add_ball_and_stick()
        if show_cell:
            if parent_basis.cell is not None:
                if all(np.max(parent_basis.cell, axis=0) > 1e-2):
                    view.add_unitcell()
        if show_axes:
            view.shape.add_arrow([-2, -2, -2], [2, -2, -2], [1, 0, 0], 0.5)
            view.shape.add_arrow([-2, -2, -2], [-2, 2, -2], [0, 1, 0], 0.5)
            view.shape.add_arrow([-2, -2, -2], [-2, -2, 2], [0, 0, 1], 0.5)
        if camera != "perspective" and camera != "orthographic":
            print("Only perspective or orthographic is permitted")
            return None
        view.camera = camera
        view.background = background
        return view

    def pos_xyz(self):
        """

        Returns:

        """
        x = self.positions[:, 0]
        y = self.positions[:, 1]
        z = self.positions[:, 2]
        return x, y, z

    def scaled_pos_xyz(self):
        """

        Returns:

        """
        xyz = self.get_scaled_positions(wrap=False)
        return xyz[:, 0], xyz[:, 1], xyz[:, 2]

    def __select_slice(self, i_dim, i_flag, dist):
        """

        Args:
            i_dim:
            i_flag:
            dist:

        Returns:

        """
        if i_dim + 1 > self.dimension:
            return True
        if i_flag == 1:
            return self.get_scaled_positions(wrap=False)[:, i_dim] < dist
        elif i_flag == 0:
            return True
        elif i_flag == -1:
            return self.get_scaled_positions(wrap=False)[:, i_dim] > 1.0 - dist

    def get_boundary_region(self, dist):
        """
        get all atoms in the boundary around the supercell which have a distance
        to the supercell boundary of less than dist

        Args:
            dist:

        Returns:

        """
        rel_coordinates = self.get_scaled_positions(wrap=False)

        dim = self.dimension
        cell = self.cell.T  # to use same definition as ASE
        a1 = cell[0]
        a2, a3 = 0, 0
        min_i, max_i = -1, 2
        iyl, iy, izl, iz = 0, 1, 0, 1
        if dim > 1:
            a2 = cell[1]
            iyl, iy = min_i, max_i
        if dim > 2:
            a3 = cell[2]
            izl, iz = min_i, max_i

        index = np.arange(len(self))
        new_coordinates = np.zeros((1, dim))
        # pbcVec = np.zeros((1, dim))
        ia_list = np.zeros((1, 1), dtype=np.int)
        for i0 in range(min_i, max_i):
            for i1 in range(iyl, iy):
                for i2 in range(izl, iz):
                    # r_vec_abs = i0 * a1 + i1 * a2 + i2 * a3
                    r_vec = np.array([i0, i1, i2][:dim])
                    select = (
                        self.__select_slice(0, i0, dist)
                        & self.__select_slice(1, i1, dist)
                        & self.__select_slice(2, i2, dist)
                    )
                    if np.linalg.norm(r_vec) > 0:
                        if len(select) > 0:
                            sel_coordinates = rel_coordinates[select] + r_vec
                            new_coordinates = np.append(
                                new_coordinates, sel_coordinates, axis=0
                            )
                            if len(sel_coordinates) > 0:
                                # rVecs = np.array(len(sel_coordinates) * [r_vec_abs])
                                # pbcVec = np.append(pbcVec, rVecs, axis=0)
                                ia_list = np.append(ia_list, index[select])
                                # print "rVec: ", i0,i1,i2,rVecs[0],index[select],select

        element_list = [self.indices[ia] for ia in ia_list[1:]]
        self._ia_bounds = ia_list[1:]
        # self._pbcVec = pbcVec[1:]
        return Atoms(
            indices=element_list,
            scaled_positions=new_coordinates[1:],
            cell=self.cell,
            dimension=len(cell),
            species=self.species,
        )

    def get_neighbors(
        self,
        num_neighbors=12,
        t_vec=True,
        include_boundary=True,
        exclude_self=True,
        tolerance=2,
        id_list=None,
        cutoff_radius=None,
        cutoff=None,
    ):
        """

        Args:
            num_neighbors (int): number of neighbors
            t_vec (bool): True: compute distance vectors
                        (pbc are automatically taken into account)
            include_boundary (bool): True: search for neighbors assuming periodic boundary conditions
                                     False is needed e.g. in plot routines to avoid showing incorrect bonds
            exclude_self (bool): include central __atom (i.e. distance = 0)
            tolerance (int): tolerance (round decimal points) used for computing neighbor shells
            id_list:
            cutoff (float/None): Upper bound of the distance to which the search must be done - by default search for
                                 upto 100 neighbors unless num_neighbors is defined explicitly.
            cutoff_radius (float/None): Upper bound of the distance to which the search must be done - by default search
                                        for upto 100 neighbors unless num_neighbors is defined explicitly.

        Returns:

            pyiron.atomistics.structure.atoms.Neighbors: Neighbors instances with the neighbor indices, distances
            and vectors

        """
        if cutoff is not None and cutoff_radius is None:
            warnings.warn(
                "Please use cutoff_radius, rather than cutoff", DeprecationWarning
            )
            cutoff_radius = cutoff
        if cutoff_radius is not None and num_neighbors == 12:
            num_neighbors = 100
        # eps = 1e-4
        i_start = 0
        if exclude_self:
            i_start = 1

        def f_ind(x):
            return x < len(self)

        num_neighbors += 1
        neighbor_obj = Neighbors()
        if not include_boundary:  # periodic boundaries are NOT included
            tree = cKDTree(self.positions)
            if cutoff_radius is None:
                neighbors = tree.query(self.positions, k=num_neighbors)
            else:
                neighbors = tree.query(
                    self.positions, k=num_neighbors, distance_upper_bound=cutoff_radius
                )

            d_lst, ind_lst, v_lst = [], [], []
            ic = 0
            for d_i, ind_i in zip(neighbors[0], neighbors[1]):
                ff = (ind_i < len(self)) & (ind_i != ic)
                ind_l = ind_i[ff]
                ind_lst.append(ind_l)
                d_lst.append(d_i[ff])
                v_lst.append(self.positions[ind_l] - self.positions[ic])
                ic += 1
            neighbor_obj.indices = ind_lst
            neighbor_obj.distances = d_lst
            neighbor_obj.vecs = v_lst
            return neighbor_obj

        # include periodic boundaries
        # translate radius in boundary layer with relative coordinates
        # TODO: introduce more rigoros definition
        radius = 3 * num_neighbors ** (1.0 / 3.0)
        rel_width = [radius / np.sqrt(np.dot(a_i, a_i)) for a_i in self.cell]
        rel_width_scalar = np.max(rel_width)

        # construct cell with additional atoms bounding original cell
        boundary_atoms = self.get_boundary_region(rel_width_scalar)
        extended_cell = self + boundary_atoms

        # build index to map boundary atoms back to original cell
        map_to_cell = np.append(np.arange(len(self)), self._ia_bounds)

        # transfer relative to absolute coordinates
        tree = cKDTree(extended_cell.positions)
        if id_list is None:
            positions = self.positions
        else:
            positions = np.array([self.positions[i] for i in id_list])
        # print ("len positions: ", len(positions))
        if cutoff_radius is None:
            neighbors = tree.query(positions, k=num_neighbors)
        else:
            neighbors = tree.query(
                positions, k=num_neighbors, distance_upper_bound=cutoff_radius
            )

        # print ("neighbors: ", neighbors)

        self.neighbor_distance = []  # neighbors[0]
        self.neighbor_distance_vec = []
        self.neighbor_index = []

        self.neighbor_shellOrder = []

        # tolerance = 2 # tolerance for round floating point

        def f_ind_ext(x):
            return x < len(extended_cell)

        neighbor_index = map(lambda x: filter(f_ind_ext, x), neighbors[1])
        num_neighbors = []
        for i, index in enumerate(neighbor_index):
            # print "i, index: ", i, index
            index = list(index)  # Filter conversion for python 3 compatibility
            nbrs_distances = neighbors[0][i][i_start : len(index)]
            # if radius:  # reduce neighborlist based on radius
            #     new_index_lst, new_dist_lst = [], []
            #     for index_red, dis_red in zip(index, nbrs_distances):
            #         if dis_red < radius:
            #             new_index_lst.append(index_red)
            #             new_dist_lst.append(dis_red)
            #     index, nbrs_distances= new_index_lst, new_dist_lst
            self.neighbor_distance.append(nbrs_distances)
            self.neighbor_index.append(map_to_cell[index][i_start:])
            u, indices = np.unique(
                np.around(nbrs_distances, decimals=tolerance), return_inverse=True
            )
            self.neighbor_shellOrder.append(
                indices + 1
            )  # this gives the shellOrder of neighboring atoms back

            if t_vec:
                nbr_dist = []
                if len(index) == 0:
                    self.neighbor_distance_vec.append(nbr_dist)
                    continue
                vec0 = self.positions[index[0]]
                for i_nbr, ind in enumerate(index[i_start:]):
                    # ind0 = map_to_cell[ind]
                    vec_r_ij = extended_cell.positions[ind] - vec0

                    dd0 = neighbors[0][i][i_nbr + i_start]
                    dd = np.sqrt(np.dot(vec_r_ij, vec_r_ij))
                    if not (dd - dd0 < 0.001):
                        raise AssertionError()
                    # if (dd - dd0 > 0.001):
                    #     print "wrong: ", vec_r_ij, dd,dd0,i_nbr,ind,ind0,i
                    #     print self.positions[ind0], extended_cell.positions[ind], vec0
                    nbr_dist.append(vec_r_ij)
                self.neighbor_distance_vec.append(nbr_dist)
            num_neighbors.append(len(index) - i_start)

        min_nbr, max_nbr = min(num_neighbors), max(num_neighbors)
        if max_nbr == num_neighbors:
            # print "neighbor distance: ", self.neighbor_distance
            raise ValueError(
                "Increase max_num_neighbors! " + str(max_nbr) + " " + str(num_neighbors)
            )
        self.min_nbr_number = min_nbr
        self.max_nbr_number = max_nbr
        neighbor_obj.distances = self.neighbor_distance
        neighbor_obj.vecs = self.neighbor_distance_vec
        neighbor_obj.indices = self.neighbor_index
        neighbor_obj.shells = self.neighbor_shellOrder
        return neighbor_obj

    def get_neighborhood(
        self,
        position,
        num_neighbors=12,
        t_vec=True,
        include_boundary=True,
        tolerance=2,
        id_list=None,
        cutoff=None,
        cutoff_radius=None,
    ):
        """

        Args:
            position: position in a box whose neighborhood information is analysed
            num_neighbors:
            t_vec (bool): True: compute distance vectors
                        (pbc are automatically taken into account)
            include_boundary (bool): True: search for neighbors assuming periodic boundary conditions
                                     False is needed e.g. in plot routines to avoid showing incorrect bonds
            tolerance (int): tolerance (round decimal points) used for computing neighbor shells
            id_list:
            cutoff (float/ None): Upper bound of the distance to which the search must be done
            cutoff_radius (float/ None): Upper bound of the distance to which the search must be done

        Returns:

            pyiron.atomistics.structure.atoms.Neighbors: Neighbors instances with the neighbor indices, distances
            and vectors

        """

        class NeighTemp(object):
            pass

        box = self.copy()
        box += box[-1]
        pos = box.positions
        pos[-1] = np.array(position)
        box.positions = pos
        neigh = box.get_neighbors(
            num_neighbors=num_neighbors,
            t_vec=t_vec,
            include_boundary=include_boundary,
            exclude_self=True,
            tolerance=tolerance,
            id_list=id_list,
            cutoff=cutoff,
            cutoff_radius=cutoff_radius,
        )
        neigh_return = NeighTemp()
        setattr(neigh_return, "distances", neigh.distances[-1])
        setattr(neigh_return, "shells", neigh.shells[-1])
        setattr(neigh_return, "vecs", neigh.vecs[-1])
        setattr(neigh_return, "indices", neigh.indices[-1])
        neigh_return.distances = neigh_return.distances[
            neigh_return.indices != len(box) - 1
        ]
        neigh_return.shells = neigh_return.shells[neigh_return.indices != len(box) - 1]
        neigh_return.vecs = np.array(neigh_return.vecs)[
            neigh_return.indices != len(box) - 1
        ]
        neigh_return.indices = neigh_return.indices[
            neigh_return.indices != len(box) - 1
        ]
        return neigh_return

    def find_neighbors_by_vector(self, vector, deviation=False, num_neighbors=96):
        """
        Args:
            vector (list/np.ndarray): vector by which positions are translated (and neighbors are searched)
            deviation (bool): whether to return distance between the expect positions and real positions
            num_neighbors (int): number of neighbors to take into account in get_neighbors

        Returns:
            np.ndarray: list of id's for the specified translation

        Example:
            a_0 = 2.832
            structure = pr.create_structure('Fe', 'bcc', a_0)
            id_list = structure.find_neighbors_by_vector([0, 0, a_0])
            # In this example, you get a list of neighbor atom id's at z+=a_0 for each atom.
            # This is particularly powerful for SSA when the magnetic structure has to be translated
            # in each direction.
        """

        neigh = self.get_neighbors(num_neighbors=num_neighbors)
        dist = np.linalg.norm(neigh.vecs-np.array(vector), axis=-1)
        if deviation:
            return neigh.indices[np.arange(len(dist)), np.argmin(dist, axis=-1)], np.min(dist, axis=-1)
        return neigh.indices[np.arange(len(dist)), np.argmin(dist, axis=-1)]

    def get_shells(self, id_list=None, max_shell=2, max_num_neighbors=100):
        """

        Args:
            id_list:
            max_shell:
            max_num_neighbors:

        Returns:

        """
        if id_list is None:
            id_list = [0]
        neighbors = self.get_neighbors(num_neighbors=max_num_neighbors, id_list=id_list)

        shells = neighbors.shells[0]
        dist = neighbors.distances[0]

        shell_dict = {}
        for i_shell in set(shells):
            if i_shell > max_shell:
                break
            shell_dict[i_shell] = np.mean(dist[shells == i_shell])
            # print ("shells: ", i_shell, shell_dict[i_shell])
        if not (max(shell_dict.keys()) == max_shell):
            raise AssertionError()
        return shell_dict

    def get_shell_matrix(
        self, shell=None, id_list=None, restraint_matrix=None, num_neighbors=100, tolerance=2
    ):
        """

        Args:
            neigh_list: user defined get_neighbors (recommended if atoms are displaced from the ideal positions)
            id_list: cf. get_neighbors
            radius: cf. get_neighbors
            num_neighbors: cf. get_neighbors
            tolerance: cf. get_neighbors
            restraint_matrix: NxN matrix with True or False, where False will remove the entries.
                              If an integer is given the sum of the chemical indices corresponding to the number will
                              be set to True and the rest to False

        Returns:
            NxN matrix with 1 for the pairs of atoms in the given shell

        """
        if shell is not None and shell<=0:
            raise ValueError("Parameter 'shell' must be an integer greater than 0")
        neigh_list = self.get_neighbors(
            num_neighbors=num_neighbors, id_list=id_list, tolerance=tolerance
        )
        Natom = len(self)
        if shell is None:
            shell_lst = np.unique(neigh_list.shells)
        else:
            shell_lst = np.array([shell]).flatten()
        if restraint_matrix is None:
            restraint_matrix = np.ones((Natom, Natom)) == 1
        elif type(restraint_matrix) == list and len(restraint_matrix) == 2:
            restraint_matrix = np.outer(
                1 * (self.get_chemical_symbols() == restraint_matrix[0]),
                1 * (self.get_chemical_symbols() == restraint_matrix[1]),
            )
            restraint_matrix = (restraint_matrix + restraint_matrix.transpose()) > 0
        shell_matrix_lst = []
        for shell in shell_lst:
            shell_matrix = np.zeros((Natom, Natom))
            for ii, ss in enumerate(neigh_list.shells):
                unique, counts = np.unique(
                    neigh_list.indices[ii][ss == np.array(shell)], return_counts=True
                )
                shell_matrix[ii][unique] = counts
            shell_matrix[restraint_matrix == False] = 0
            shell_matrix_lst.append(shell_matrix)
        if len(shell_matrix_lst)==1:
            return shell_matrix_lst[0]
        else:
            return shell_matrix_lst

    def get_shell_radius(self, shell=1, id_list=None):
        """

        Args:
            shell:
            id_list:

        Returns:

        """
        if id_list is None:
            id_list = [0]
        shells = self.get_shells(id_list=id_list, max_shell=shell + 1)
        return np.mean(list(shells.values())[shell - 1 :])

    def occupy_lattice(self, **qwargs):
        """
        Replaces specified indices with a given species
        """
        new_species = list(np.array(self.species).copy())
        new_indices = np.array(self.indices.copy())
        for key, i_list in qwargs.items():
            el = self._pse.element(key)
            if el.Abbreviation not in [spec.Abbreviation for spec in new_species]:
                new_species.append(el)
                new_indices[i_list] = len(new_species) - 1
            else:
                index = np.argwhere(np.array(new_species) == el).flatten()
                new_indices[i_list] = index
        delete_species_indices = list()
        retain_species_indices = list()
        for i, el in enumerate(new_species):
            if len(np.argwhere(new_indices == i).flatten()) == 0:
                delete_species_indices.append(i)
            else:
                retain_species_indices.append(i)
        for i in delete_species_indices:
            new_indices[new_indices >= i] += -1
        new_species = np.array(new_species)[retain_species_indices]
        self.set_species(new_species)
        self.indices = new_indices

    def cluster_analysis(
        self, id_list, neighbors=None, radius=None, return_cluster_sizes=False
    ):
        """

        Args:
            id_list:
            neighbors:
            radius:
            return_cluster_sizes:

        Returns:

        """
        if neighbors is None:
            if radius is None:
                radius = self.get_shell_radius()
                # print "radius: ", radius
            neighbors = self.get_neighbors(radius, t_vec=False)
        self._neighbor_index = neighbors.indices
        self._cluster = [0] * len(self)
        c_count = 1
        # element_list = self.get_atomic_numbers()
        for ia in id_list:
            # el0 = element_list[ia]
            nbrs = self._neighbor_index[ia]
            # print ("nbrs: ", ia, nbrs)
            if self._cluster[ia] == 0:
                self._cluster[ia] = c_count
                self.__probe_cluster(c_count, nbrs, id_list)
                c_count += 1

        cluster = np.array(self._cluster)
        cluster_dict = {
            i_c: np.where(cluster == i_c)[0].tolist() for i_c in range(1, c_count)
        }
        if return_cluster_sizes:
            sizes = [self._cluster.count(i_c + 1) for i_c in range(c_count - 1)]
            return cluster_dict, sizes

        return cluster_dict  # sizes

    def __probe_cluster(self, c_count, neighbors, id_list):
        """

        Args:
            c_count:
            neighbors:
            id_list:

        Returns:

        """
        for nbr_id in neighbors:
            if self._cluster[nbr_id] == 0:
                if nbr_id in id_list:  # TODO: check also for ordered structures
                    self._cluster[nbr_id] = c_count
                    nbrs = self._neighbor_index[nbr_id]
                    self.__probe_cluster(c_count, nbrs, id_list)

    # TODO: combine with corresponding routine in plot3d
    def get_bonds(self, radius=None, max_shells=None, prec=0.1, num_neighbors=20):
        """

        Args:
            radius:
            max_shells:
            prec: minimum distance between any two clusters (if smaller considered to be single cluster)
            num_neighbors:

        Returns:

        """

        def get_cluster(dist_vec, ind_vec, prec=prec):
            ind_where = np.where(np.diff(dist_vec) > prec)[0] + 1
            ind_vec_cl = [np.sort(group) for group in np.split(ind_vec, ind_where)]
            dist_vec_cl = [np.mean(group) for group in np.split(dist_vec, ind_where)]
            return ind_vec_cl, dist_vec_cl

        neighbors = self.get_neighbors(
            cutoff_radius=radius, num_neighbors=num_neighbors
        )

        dist = neighbors.distances
        ind = neighbors.indices
        el_list = self.get_chemical_symbols()

        ind_shell = []
        for i_a, (d, i) in enumerate(zip(dist, ind)):
            id_list, dist_lst = get_cluster(d[d < radius], i[d < radius])
            # print ("id: ", d[d<radius], id_list, dist_lst)
            ia_shells_dict = {}
            for i_shell_list in id_list:
                ia_shell_dict = {}
                for i_s in i_shell_list:
                    el = el_list[i_s]
                    if el not in ia_shell_dict:
                        ia_shell_dict[el] = []
                    ia_shell_dict[el].append(i_s)
                for el, ia_lst in ia_shell_dict.items():
                    if el not in ia_shells_dict:
                        ia_shells_dict[el] = []
                    if max_shells is not None:
                        if len(ia_shells_dict[el]) + 1 > max_shells:
                            continue
                    ia_shells_dict[el].append(ia_lst)
            ind_shell.append(ia_shells_dict)
        return ind_shell

    # spglib calls
    def get_symmetry(
        self, use_magmoms=False, use_elements=True, symprec=1e-5, angle_tolerance=-1.0
    ):
        """

        Args:
            use_magmoms:
            use_elements: True or False. If False, chemical elements will be ignored
            symprec:
            angle_tolerance:

        Returns:


        """
        lattice = np.array(self.get_cell().T, dtype="double", order="C")
        positions = np.array(
            self.get_scaled_positions(wrap=False), dtype="double", order="C"
        )
        if use_elements:
            numbers = np.array(self.get_atomic_numbers(), dtype="intc")
        else:
            numbers = np.ones_like(self.get_atomic_numbers(), dtype="intc")
        if use_magmoms:
            magmoms = self.get_initial_magnetic_moments()
            return spglib.get_symmetry(
                cell=(lattice, positions, numbers, magmoms),
                symprec=symprec,
                angle_tolerance=angle_tolerance,
            )
        else:
            return spglib.get_symmetry(
                cell=(lattice, positions, numbers),
                symprec=symprec,
                angle_tolerance=angle_tolerance,
            )

    def symmetrize_vectors(
        self, vectors, force_update=False, use_magmoms=False, use_elements=True, symprec=1e-5, angle_tolerance=-1.0
    ):
        """
        Symmetrization of natom x 3 vectors according to box symmetries

        Args:
            vectors (ndarray/list): natom x 3 array to symmetrize
            force_update (bool): whether to update the symmetry info
            use_magmoms (bool): cf. get_symmetry
            use_elements (bool): cf. get_symmetry
            symprec (float): cf. get_symmetry
            angle_tolerance (float): cf. get_symmetry

        Returns:
            (np.ndarray) symmetrized vectors
        """
        vectors = np.array(vectors).reshape(-1, 3)
        if vectors.shape != self.positions.shape:
            print(vectors.shape, self.positions.shape)
            raise ValueError('Vector must be a natom x 3 array: {} != {}'.format(vectors.shape, self.positions.shape))
        if self._symmetry_dataset is None or force_update:
            symmetry = self.get_symmetry(use_magmoms=use_magmoms, use_elements=use_elements,
                                         symprec=symprec, angle_tolerance=angle_tolerance)
            scaled_positions = self.get_scaled_positions(wrap=False)
            symmetry['indices'] = []
            for rot,tra in zip(symmetry['rotations'], symmetry['translations']):
                positions = np.einsum('ij,nj->ni', rot, scaled_positions)+tra
                positions -= np.floor(positions+1.0e-2)
                vec = np.where(np.linalg.norm(positions[np.newaxis, :, :]-scaled_positions[:, np.newaxis, :], axis=-1)<=1.0e-4)
                symmetry['indices'].append(vec[1])
            symmetry['indices'] = np.array(symmetry['indices'])
            self._symmetry_dataset = symmetry
        return np.einsum('ijk,ink->nj', self._symmetry_dataset['rotations'],
                         vectors[self._symmetry_dataset['indices']])/len(self._symmetry_dataset['rotations'])

    def group_points_by_symmetry(self, points):
        """
            This function classifies the points into groups according to the box symmetry given by spglib.

            Args:
                points: (np.array/list) nx3 array which contains positions

            Returns: list of arrays containing geometrically equivalent positions

            It is possible that the original points are not found in the returned list, as the positions outsie
            the box will be projected back to the box.
        """
        struct_copy = self.copy()
        points = np.array(points).reshape(-1, 3)
        struct_copy += Atoms(elements=len(points) * ["Hs"], positions=points)
        struct_copy.center_coordinates_in_unit_cell()
        group_IDs = struct_copy.get_symmetry()["equivalent_atoms"][
            struct_copy.select_index("Hs")
        ]
        return [
            np.round(points[group_IDs == ID], decimals=8) for ID in np.unique(group_IDs)
        ]

    def _get_voronoi_vertices(self, minimum_dist=0.1):
        """
            This function gives the positions of Voronoi vertices
            This function does not work if there are Hs atoms in the box

            Args:
                minimum_dist: Minimum distance between two Voronoi vertices to be considered as one

            Returns: Positions of Voronoi vertices, box

        """
        vor = Voronoi(
            self.repeat(3 * [2]).positions
        )  # Voronoi package does not have periodic boundary conditions
        b_cell_inv = np.linalg.inv(self.cell)
        voro_vert = vor.vertices
        for ind, v in enumerate(voro_vert):
            pos = np.mean(
                voro_vert[(np.linalg.norm(voro_vert - v, axis=-1) < minimum_dist)],
                axis=0,
            )  # Find all points which are within minimum_dist
            voro_vert[(np.linalg.norm(voro_vert - v, axis=-1) < 0.5)] = np.array(
                3 * [-10]
            )  # Mark atoms to be deleted afterwards
            voro_vert[ind] = pos
        voro_vert = voro_vert[np.min(voro_vert, axis=-1) > -5]

        voro_vert = np.dot(b_cell_inv.T, voro_vert.T).T  # get scaled positions
        voro_vert = voro_vert[
            (np.min(voro_vert, axis=-1) > 0.499) & (np.max(voro_vert, axis=-1) < 1.501)
        ]
        voro_vert = np.dot(self.cell.T, voro_vert.T).T  # get true positions

        box_copy = self.copy()
        new_atoms = Atoms(cell=self.cell, symbols=["Hs"]).repeat([len(voro_vert), 1, 1])
        box_copy += new_atoms

        pos_total = np.append(self.positions, voro_vert)
        pos_total = pos_total.reshape(-1, 3)
        box_copy.positions = pos_total

        box_copy.center_coordinates_in_unit_cell()

        neigh = (
            box_copy.get_neighbors()
        )  # delete all atoms which lie within minimum_dist (including periodic boundary conditions)
        while (
            len(
                np.array(neigh.indices).flatten()[
                    np.array(neigh.distances).flatten() < minimum_dist
                ]
            )
            != 0
        ):
            del box_copy[
                np.array(neigh.indices).flatten()[
                    np.array(neigh.distances).flatten() < minimum_dist
                ][0]
            ]
            neigh = box_copy.get_neighbors()
        return pos_total, box_copy

    def get_equivalent_voronoi_vertices(
        self, return_box=False, minimum_dist=0.1, symprec=1e-5, angle_tolerance=-1.0
    ):
        """
            This function gives the positions of spatially equivalent Voronoi vertices in lists, which
            most likely represent interstitial points or vacancies (along with other high symmetry points)
            Each list item contains an array of positions which are spacially equivalent.
            This function does not work if there are Hs atoms in the box

            Args:
                return_box: True, if the box containing atoms on the positions of Voronoi vertices
                            should be returned (which are represented by Hs atoms)
                minimum_dist: Minimum distance between two Voronoi vertices to be considered as one

            Returns: List of numpy array positions of spacially equivalent Voronoi vertices

        """

        _, box_copy = self._get_voronoi_vertices(minimum_dist=minimum_dist)
        list_positions = []
        sym = box_copy.get_symmetry(symprec=symprec, angle_tolerance=angle_tolerance)
        for ind in set(sym["equivalent_atoms"][box_copy.select_index("Hs")]):
            list_positions.append(box_copy.positions[sym["equivalent_atoms"] == ind])
        if return_box:
            return list_positions, box_copy
        else:
            return list_positions

    def get_equivalent_points(self, points, use_magmoms=False, use_elements=True, symprec=1e-5, angle_tolerance=-1.0):
        """

        Args:
            points (list/ndarray): 3d vector
            use_magmoms (bool): cf. get_symmetry()
            use_elements (bool): cf. get_symmetry()
            symprec (float): cf. get_symmetry()
            angle_tolerance (float): cf. get_symmetry()

        Returns:
            (ndarray): array of equivalent points with respect to box symmetries
        """
        symmetry_operations = self.get_symmetry(use_magmoms=use_magmoms,
                                                use_elements=use_elements,
                                                symprec=symprec,
                                                angle_tolerance=angle_tolerance)
        R = symmetry_operations['rotations']
        t = symmetry_operations['translations']
        x = np.einsum('jk,j->k', np.linalg.inv(self.cell), points)
        x = np.einsum('nxy,y->nx', R, x)+t
        x -= np.floor(x)
        dist = x[:,np.newaxis]-x[np.newaxis,:]
        w, v = np.where(np.linalg.norm(dist-np.rint(dist), axis=-1)<symprec)
        x = np.delete(x, w[v<w], axis=0)
        x = np.einsum('ji,mj->mi', self.cell, x)
        return x

    def get_symmetry_dataset(self, symprec=1e-5, angle_tolerance=-1.0):
        """

        Args:
            symprec:
            angle_tolerance:

        Returns:

        https://atztogo.github.io/spglib/python-spglib.html
        """
        lattice = np.array(self.get_cell().T, dtype="double", order="C")
        positions = np.array(
            self.get_scaled_positions(wrap=False), dtype="double", order="C"
        )
        numbers = np.array(self.get_atomic_numbers(), dtype="intc")
        return spglib.get_symmetry_dataset(
            cell=(lattice, positions, numbers),
            symprec=symprec,
            angle_tolerance=angle_tolerance,
        )

    def get_spacegroup(self, symprec=1e-5, angle_tolerance=-1.0):
        """

        Args:
            symprec:
            angle_tolerance:

        Returns:

        https://atztogo.github.io/spglib/python-spglib.html
        """
        lattice = np.array(self.get_cell(), dtype="double", order="C")
        positions = np.array(
            self.get_scaled_positions(wrap=False), dtype="double", order="C"
        )
        numbers = np.array(self.get_atomic_numbers(), dtype="intc")
        space_group = spglib.get_spacegroup(
            cell=(lattice, positions, numbers),
            symprec=symprec,
            angle_tolerance=angle_tolerance,
        ).split()
        if len(space_group) == 1:
            return {"Number": ast.literal_eval(space_group[0])}
        else:
            return {
                "InternationalTableSymbol": space_group[0],
                "Number": ast.literal_eval(space_group[1]),
            }

    def refine_cell(self, symprec=1e-5, angle_tolerance=-1.0):
        """

        Args:
            symprec:
            angle_tolerance:

        Returns:

        https://atztogo.github.io/spglib/python-spglib.html
        """
        lattice = np.array(self.get_cell().T, dtype="double", order="C")
        positions = np.array(
            self.get_scaled_positions(wrap=False), dtype="double", order="C"
        )
        numbers = np.array(self.get_atomic_numbers(), dtype="intc")
        cell, coords, el = spglib.refine_cell(
            cell=(lattice, positions, numbers),
            symprec=symprec,
            angle_tolerance=angle_tolerance,
        )

        return Atoms(
            symbols=list(self.get_chemical_symbols()), positions=coords, cell=cell
        )

    def get_primitive_cell(self, symprec=1e-5, angle_tolerance=-1.0):
        """

        Args:
            symprec:
            angle_tolerance:

        Returns:

        """
        el_dict = {}
        for el in set(self.get_chemical_elements()):
            el_dict[el.AtomicNumber] = el
        lattice = np.array(self.get_cell().T, dtype="double", order="C")
        positions = np.array(
            self.get_scaled_positions(wrap=False), dtype="double", order="C"
        )
        numbers = np.array(self.get_atomic_numbers(), dtype="intc")
        cell, coords, atomic_numbers = spglib.find_primitive(
            cell=(lattice, positions, numbers),
            symprec=symprec,
            angle_tolerance=angle_tolerance,
        )
        # print atomic_numbers, type(atomic_numbers)
        el_lst = [el_dict[i_a] for i_a in atomic_numbers]

        # convert lattice vectors to standard (experimental feature!) TODO:
        red_structure = Atoms(elements=el_lst, scaled_positions=coords, cell=cell)
        space_group = red_structure.get_spacegroup(symprec)["Number"]
        # print "space group: ", space_group
        if space_group == 225:  # fcc
            alat = np.max(cell[0])
            amat_fcc = alat * np.array([[1, 0, 1], [1, 1, 0], [0, 1, 1]])

            red_structure.cell = amat_fcc
        return red_structure

    def get_ir_reciprocal_mesh(
        self,
        mesh,
        is_shift=np.zeros(3, dtype="intc"),
        is_time_reversal=True,
        symprec=1e-5,
    ):
        """

        Args:
            mesh:
            is_shift:
            is_time_reversal:
            symprec:

        Returns:

        """
        mapping, mesh_points = spglib.get_ir_reciprocal_mesh(
            mesh=mesh,
            cell=self,
            is_shift=is_shift,
            is_time_reversal=is_time_reversal,
            symprec=symprec,
        )
        return mapping, mesh_points

    def get_majority_species(self, return_count=False):
        """
        This function returns the majority species and their number in the box

        Returns:
            number of atoms of the majority species, chemical symbol and chemical index

        """
        el_dict = self.get_number_species_atoms()
        el_num = list(el_dict.values())
        el_name = list(el_dict.keys())
        if np.sum(np.array(el_num) == np.max(el_num)) > 1:
            warnings.warn("There are more than one majority species")
        symbol_to_index = dict(
            zip(self.get_chemical_symbols(), self.get_chemical_indices())
        )
        max_index = np.argmax(el_num)
        return {
            "symbol": el_name[max_index],
            "count": int(np.max(el_num)),
            "index": symbol_to_index[el_name[max_index]],
        }

    def close(self):
        # TODO: implement
        pass

    def get_voronoi_volume(self):
        """

        Returns:

        """
        warnings.warn(
            "This function doesn't account for periodic boundary conditions. Call "
            "`analyse_ovito_voronoi_volume` instead. This is what will now be returned.",
            DeprecationWarning,
        )
        return self.analyse_ovito_voronoi_volume()

    def extend(self, other):
        """
        Extend atoms object by appending atoms from *other*. (Extending the ASE function)

        Args:
            other (pyiron.atomistics.structure.atoms.Atoms): Structure to append

        Returns:
            pyiron.atomistics.structure.atoms.Atoms: The extended structure

        """
        old_indices = self.indices
        if isinstance(other, Atom):
            other = self.__class__([other])
        new_indices = other.indices
        super(Atoms, self).extend(other=other)
        if isinstance(other, Atoms):
            sum_atoms = self
            # sum_atoms = copy(self)
            sum_atoms._tag_list = sum_atoms._tag_list + other._tag_list
            sum_atoms.indices = np.append(sum_atoms.indices, other.indices)
            new_species_lst = copy(sum_atoms.species)
            ind_conv = {}
            for ind_old, el in enumerate(other.species):
                if el.Abbreviation in sum_atoms._store_elements.keys():
                    ind_new = sum_atoms._species_to_index_dict[
                        sum_atoms._store_elements[el.Abbreviation]
                    ]
                    ind_conv[ind_old] = ind_new
                else:
                    new_species_lst.append(el)
                    sum_atoms._store_elements[el.Abbreviation] = el
                    ind_conv[ind_old] = len(new_species_lst) - 1

            for key, val in ind_conv.items():
                new_indices[new_indices == key] = val + 1000
            new_indices = np.mod(new_indices, 1000)
            sum_atoms.indices[len(old_indices):] = new_indices
            sum_atoms.set_species(new_species_lst)
            if not len(set(sum_atoms.indices)) == len(sum_atoms.species):
                raise ValueError("Adding the atom instances went wrong!")
        return self

    __iadd__ = extend

    def __copy__(self):
        """
        Copies the atoms object

        Returns:
            atoms_new: A copy of the object

        """
        atoms_new = Atoms()
        for key, val in self.__dict__.items():
            if key not in ["_pse"]:
                # print ('copy: ', key)
                atoms_new.__dict__[key] = copy(val)

        return atoms_new

    def __delitem__(self, key):
        if isinstance(key, (int, np.integer)):
            key = [key]
        key = np.array(key).flatten()
        new_length = len(self) - len(key)
        super(Atoms, self).__delitem__(key)
        self.indices = np.delete(self.indices, key, axis=0)
        del self._tag_list[key]
        self._tag_list._length = new_length
        deleted_species_indices = list()
        retain_species_indices = list()
        new_indices = self.indices.copy()
        for i, el in enumerate(self.species):
            if len(self.select_index(el)) == 0:
                deleted_species_indices.append(i)
                new_indices[new_indices >= i] += -1
            else:
                retain_species_indices.append(i)
        new_species = np.array(self.species).copy()[retain_species_indices]
        self.set_species(new_species)
        self.indices = new_indices

    def __eq__(self, other):
        if not (isinstance(other, Atoms)):
            raise AssertionError()
        conditions = []
        for a_1, a_2 in zip(self, other):
            conditions.append(a_1 == a_2)
        conditions.append(np.alltrue(self.pbc == other.pbc))
        return all(conditions)

    def __ne__(self, other):
        return not self == other

    def __getitem__(self, item):
        new_dict = dict()
        if isinstance(item, int):
            for key, value in self._tag_list.items():
                if item < len(value):
                    if value[item] is not None:
                        new_dict[key] = value[item]
            element = self.species[self.indices[item]]
            index = item
            position = self.positions[item]
            return Atom(
                element=element,
                position=position,
                pse=self._pse,
                index=index,
                atoms=self,
                **new_dict
            )

        new_array = super(Atoms, self).__getitem__(item)
        new_indices = self.indices[item].copy()
        new_species_indices, new_proper_indices = np.unique(
            new_indices, return_inverse=True
        )
        new_species = [self.species[ind] for ind in new_species_indices]
        new_array.set_species(new_species)
        new_array.indices = new_proper_indices
        new_array._tag_list = self._tag_list[item]
        # new_array._tag_list._length = self._tag_list._length
        new_array._tag_list._length = len(new_array)
        if isinstance(new_array, Atom):
            natoms = len(self)
            if item < -natoms or item >= natoms:
                raise IndexError("Index out of range.")
            new_array.index = item
        return new_array

    def __getattr__(self, item):
        if item in self._tag_list.keys():
            return self._tag_list._lists[item]
        return object.__getattribute__(self, item)

    # def __len__(self):
    #     return len(self.indices)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        if len(self) == 0:
            return "[]"
        out_str = ""
        for el, pos in zip(self.get_chemical_symbols(), self.positions):
            out_str += el + ": " + str(pos) + "\n"
        if len(self.get_tags()) > 0:
            tags = self.get_tags()
            out_str += "tags: \n"  # + ", ".join(tags) + "\n"
            for tag in tags:
                out_str += (
                    "    " + str(tag) + ": " + self._tag_list[tag].__str__() + "\n"
                )
        if self.cell is not None:
            out_str += "pbc: " + str(self.pbc) + "\n"
            out_str += "cell: \n"
            out_str += str(self.cell) + "\n"
        return out_str

    def __setitem__(self, key, value):
        if isinstance(key, (int, np.integer)):
            old_el = self.species[self.indices[key]]
            if isinstance(value, (str, np.str, np.str_)):
                el = PeriodicTable().element(value)
            elif isinstance(value, ChemicalElement):
                el = value
            else:
                raise TypeError("value should either be a string or a ChemicalElement.")
            if el != old_el:
                new_species = np.array(self.species).copy()
                if len(self.select_index(old_el)) == 1:
                    if el.Abbreviation not in [
                        spec.Abbreviation for spec in new_species
                    ]:
                        new_species[self.indices[key]] = el
                        self.set_species(list(new_species))
                    else:
                        el_list = np.array([sp.Abbreviation for sp in new_species])
                        ind = np.argwhere(el_list == el.Abbreviation).flatten()[-1]
                        remove_index = self.indices[key]
                        new_species = list(new_species)
                        del new_species[remove_index]
                        self.indices[key] = ind
                        self.indices[self.indices > remove_index] -= 1
                        self.set_species(new_species)
                else:
                    if el.Abbreviation not in [
                        spec.Abbreviation for spec in new_species
                    ]:
                        new_species = list(new_species)
                        new_species.append(el)
                        self.set_species(new_species)
                        self.indices[key] = len(new_species) - 1
                    else:
                        el_list = np.array([sp.Abbreviation for sp in new_species])
                        ind = np.argwhere(el_list == el.Abbreviation).flatten()[-1]
                        self.indices[key] = ind
        elif isinstance(key, slice) or isinstance(key, (list, tuple, np.ndarray)):
            if not isinstance(key, slice):
                if hasattr(key, "__len__"):
                    if len(key) == 0:
                        return
            else:
                # Generating the correct numpy array from a slice input
                if key.start is None:
                    start_val = 0
                elif key.start < 0:
                    start_val = key.start + len(self)
                else:
                    start_val = key.start
                if key.stop is None:
                    stop_val = len(self)
                elif key.stop < 0:
                    stop_val = key.stop + len(self)
                else:
                    stop_val = key.stop
                if key.step is None:
                    step_val = 1
                else:
                    step_val = key.step
                key = np.arange(start_val, stop_val, step_val)
            if isinstance(value, (str, np.str, np.str_, int, np.integer)):
                el = PeriodicTable().element(value)
            elif isinstance(value, ChemicalElement):
                el = value
            else:
                raise ValueError(
                    "The value assigned should be a string, integer or a ChemicalElement instance"
                )
            replace_list = list()
            new_species = list(np.array(self.species).copy())

            for sp in self.species:
                replace_list.append(
                    np.array_equal(
                        np.sort(self.select_index(sp)),
                        np.sort(np.intersect1d(self.select_index(sp), key)),
                    )
                )
            if el.Abbreviation not in [spec.Abbreviation for spec in new_species]:
                if not any(replace_list):
                    new_species.append(el)
                    self.set_species(new_species)
                    self.indices[key] = len(new_species) - 1
                else:
                    replace_ind = np.where(replace_list)[0][0]
                    new_species[replace_ind] = el
                    if len(np.where(replace_list)[0]) > 1:
                        for ind in replace_list[1:]:
                            del new_species[ind]
                    self.set_species(new_species)
                    self.indices[key] = replace_ind
            else:
                el_list = np.array([sp.Abbreviation for sp in new_species])
                ind = np.argwhere(el_list == el.Abbreviation).flatten()[-1]
                if not any(replace_list):
                    self.set_species(new_species)
                    self.indices[key] = ind
                else:
                    self.indices[key] = ind
                    delete_indices = list()
                    new_indices = self.indices.copy()
                    for i, rep in enumerate(replace_list):
                        if i != ind and rep:
                            delete_indices.append(i)
                            # del new_species[i]
                            new_indices[new_indices >= i] -= 1
                    self.indices = new_indices.copy()
                    new_species = np.array(new_species)[
                        np.setdiff1d(np.arange(len(new_species)), delete_indices)
                    ].tolist()
                    self.set_species(new_species)
        else:
            raise NotImplementedError()

    __mul__ = repeat

    def __imul__(self, vec):
        if isinstance(vec, int):
            vec = [vec] * self.dimension

        if not (len(vec) == self.dimension):
            raise ValueError("The repeat vector and the dimensions don't match")
        initial_length = len(self)
        i_vec = np.array([vec[0], 1, 1])
        if self.dimension > 1:
            i_vec[1] = vec[1]
        if self.dimension > 2:
            i_vec[2] = vec[2]

        if not self.dimension == 3:
            raise NotImplementedError()
        mx, my, mz = i_vec
        # Our position repeat algorithm is faster than ASE (no nested loops)
        nx_lst, ny_lst, nz_lst = np.arange(mx), np.arange(my), np.arange(mz)
        positions = self.get_scaled_positions(wrap=False)
        lat = np.array(np.meshgrid(nx_lst, ny_lst, nz_lst)).T.reshape(-1, 3)
        lat_new = np.repeat(lat, len(positions), axis=0)
        new_positions = np.tile(positions, (len(lat), 1)) + lat_new
        new_positions /= np.array(i_vec)
        self.set_cell((self.cell.T * np.array(vec)).T, scale_atoms=True)
        # ASE compatibility
        for name, a in self.arrays.items():
            self.arrays[name] = np.tile(a, (np.product(vec),) + (1, ) * (len(a.shape) - 1))
        self.arrays["positions"] = np.dot(new_positions, self.cell)
        self.indices = np.tile(self.indices, len(lat))
        self._tag_list._length = len(self)
        scale = i_vec[0] * i_vec[1] * i_vec[2]
        for tag in self._tag_list.keys():
            self._tag_list[tag] *= scale
        # Repeating ASE constraints
        if self.constraints is not None:
            self.constraints = [c.repeat(vec, initial_length) for c in self.constraints]
        return self

    @staticmethod
    def convert_formula(elements):
        """

        Args:
            elements:

        Returns:

        """
        el_list = []
        num_list = ""
        for i, char in enumerate(elements):
            is_last = i == len(elements) - 1
            if len(num_list) > 0:
                if (not char.isdigit()) or is_last:
                    el_fac = ast.literal_eval(num_list) * el_list[-1]
                    for el in el_fac[1:]:
                        el_list.append(el)
                    num_list = ""

            if char.isupper():
                el_list.append(char)
            elif char.islower():
                el_list[-1] += char
            elif char.isdigit():
                num_list += char

            if len(num_list) > 0:
                # print "num_list: ", el_list, num_list, el_list[-1], (not char.isdigit()) or is_last
                if (not char.isdigit()) or is_last:
                    el_fac = ast.literal_eval(num_list) * [el_list[-1]]
                    # print "el_fac: ", el_fac
                    for el in el_fac[1:]:
                        el_list.append(el)
                    num_list = ""

        return el_list

    # ASE compatibility
    def write(self, filename, format=None, **kwargs):
        """
        Write atoms object to a recognized file format using ase parsers.

        see ase.io.write for formats.
        kwargs are passed to ase.io.write.
        """
        pyiron_to_ase(self).write(filename=filename, format=format, **kwargs)
    
    @staticmethod
    def get_calculator():
        return None

    # def get_cell(self, complete=False):
    #     """Get the three unit cell vectors as a 3x3 ndarray."""
    #     if complete:
    #         return complete_cell(self._cell)
    #     else:
    #         return self._cell.copy()

    def get_distance(self, a0, a1, mic=True, vector=False):
        """
        Return distance between two atoms.

        Use mic=True to use the Minimum Image Convention.
        vector=True gives the distance vector (from a0 to a1).

        Args:
            a0: position or atom ID
            a1: position or atom ID
            mic: minimum image convention (True if periodic boundary conditions should be considered)
            vector: True, if instead of distnce the vector connecting the two positions should be returned

        Returns: distance or vectors in length unit

        """
        from ase.geometry import find_mic

        positions = self.positions
        if isinstance(a0, list) or isinstance(a0, np.ndarray):
            if not (len(a0) == 3):
                raise AssertionError()
            a0 = np.array(a0)
        else:
            a0 = positions[a0]
        if isinstance(a1, list) or isinstance(a1, np.ndarray):
            if not (len(a1) == 3):
                raise AssertionError()
            a1 = np.array(a1)
        else:
            a1 = positions[a1]
        distance = np.array([a1 - a0])
        if mic:
            distance, d_len = find_mic(distance, self.cell, self.pbc)
        else:
            d_len = np.array([np.sqrt((distance ** 2).sum())])
        if vector:
            return distance[0]

        return d_len[0]

    def get_distances(self, a0=None, a1=None, mic=True, vector=False):
        """
        Return distance matrix of every position in p1 with every position in p2

        Args:
            a0 (numpy.ndarray/list): Nx3 array of positions
            a1 (numpy.ndarray/list): Nx3 array of positions
            mic (bool): minimum image convention
            vector (bool): return vectors instead of distances

        Returns:
            numpy.ndarray NxN if vector=False and NxNx3 if vector=True

        if a1 is not set, it is assumed that distances between all positions in a0 are desired. a1 will be set to a0 in this case.
        if both a0 and a1 are not set, the distances between all atoms in the box are returned

        Use mic to use the minimum image convention.

        Learn more about get_distances from the ase website:
        https://wiki.fysik.dtu.dk/ase/ase/geometry.html#ase.geometry.get_distances
        """
        if a0 is None and a1 is not None:
            a0 = a1
            a1 = None
        if a0 is None:
            a0 = self.positions
        a0 = np.array(a0).reshape(-1, 3)
        if a1 is not None:
            a1 = np.array(a1).reshape(-1, 3)
        if mic:
            vec, dist = get_distances(a0, a1, cell=self.cell, pbc=self.pbc)
        else:
            vec, dist = get_distances(a0, a1)
        if vector:
            return vec
        else:
            return dist

    def get_distance_matrix(self, mic=True, vector=False):
        """
        Return distances between all atoms in a matrix. cf. get_distance
        """
        warnings.warn(
            "get_distance_matrix is deprecated. Use get_distances instead",
            DeprecationWarning,
        )
        return self.get_distances(mic=mic, vector=vector)

    def get_constraint(self):
        if "selective_dynamics" in self._tag_list._lists.keys():
            from ase.constraints import FixAtoms

            return FixAtoms(
                indices=np.array(
                    [
                        atom_ind
                        for atom_ind in range(len(self))
                        if any(self.selective_dynamics[atom_ind])
                    ]
                )
            )
        else:
            return None

    def set_constraint(self, constraint=None):
        super(Atoms, self).set_constraint(constraint)
        if constraint is not None:
            if constraint.todict()["name"] != "FixAtoms":
                raise ValueError("Only FixAtoms is supported as ASE compatible constraint.")
            if "selective_dynamics" not in self._tag_list._lists.keys():
                self.add_tag(selective_dynamics=None)
            for atom_ind in range(len(self)):
                if atom_ind in constraint.index:
                    self.selective_dynamics[atom_ind] = [True, True, True]
                else:
                    self.selective_dynamics[atom_ind] = [False, False, False]

    def apply_strain(self, epsilon, return_box=False):
        """
            Args:
                epsilon (float/list/ndarray): epsilon matrix. If a single number is set, the same strain
                                              is applied in each direction. If a 3-dim vector is set, it
                                              will be multiplied by a unit matrix.
                return_box (bool): whether to return a box. If set to True, only the returned box will 
                                   have the desired strain and the original box will stay unchanged.
        """
        epsilon = np.array([epsilon]).flatten()
        if len(epsilon) == 3 or len(epsilon) == 1:
            epsilon = epsilon*np.eye(3)
        epsilon.reshape(3,3)
        if epsilon.min() < -1.0:
            raise ValueError("Strain value too negative")
        if return_box:
            structure_copy = self.copy()
        else:
            structure_copy = self
        cell = structure_copy.cell.copy()
        cell = np.matmul(epsilon+np.eye(3), cell)
        structure_copy.set_cell(cell, scale_atoms=True)
        if return_box:
            return structure_copy

    def get_spherical_coordinates(self, x=None):
        """
            Args:
                x (list/ndarray): coordinates to transform. If not set, the positions
                                  in structure will be returned.

            Returns:
                array in spherical coordinates
        """
        if x is None:
            x = self.positions.copy()
        x = np.array(x).reshape(-1, 3)
        r = np.linalg.norm(x, axis=-1)
        phi = np.arctan2(x[:,2], x[:,1])
        theta = np.arctan2(np.linalg.norm(x[:,:2], axis=-1), x[:,2])
        return np.stack((r, theta, phi), axis=-1)

    def get_initial_magnetic_moments(self):
        """
        Get array of initial magnetic moments.

        Returns:
            numpy.array()
        """
        if "spin" in self._tag_list._lists.keys():
            return np.asarray(self.spin.list())
        else:
            spin_lst = [
                element.tags["spin"] if "spin" in element.tags.keys() else None
                for element in self.get_chemical_elements()
            ]
            if any(spin_lst):
                if (
                    isinstance(spin_lst, str)
                    or (
                        isinstance(spin_lst, (list, np.ndarray))
                        and isinstance(spin_lst[0], str)
                    )
                ) and "[" in list(set(spin_lst))[0]:
                    return np.array(
                        [
                            [
                                float(spin_dir)
                                for spin_dir in spin.replace("[", "")
                                .replace("]", "")
                                .replace(",", "")
                                .split()
                            ]
                            if spin
                            else [0.0, 0.0, 0.0]
                            for spin in spin_lst
                        ]
                    )
                elif isinstance(spin_lst, (list, np.ndarray)):
                    return np.array(spin_lst)
                else:
                    return np.array([float(spin) if spin else 0.0 for spin in spin_lst])
            else:
                return np.array([None] * len(self))

    def set_initial_magnetic_moments(self, magmoms):
        """
        Set array of initial magnetic moments.

        Args:
            magmoms (numpy.array()):
        """
        if magmoms is not None:
            if len(magmoms) != len(self):
                raise ValueError("magmons can be collinear or non-collinear.")
            for ind, element in enumerate(self.get_chemical_elements()):
                if "spin" in element.tags.keys():
                    self[ind] = element.Parent
            if "spin" not in self._tag_list._lists.keys():
                self.add_tag(spin=None)
            for ind, spin in enumerate(magmoms):
                self.spin[ind] = spin

    def pop(self, i=-1):
        """
        Remove and return atom at index *i* (default last).

        Args:
            i:

        Returns:

        """
        atom = self[i]
        atom.cut_reference_to_atoms()
        del self[i]
        return atom

    def rotate(
        self, vector, angle=None, center=(0, 0, 0), rotate_cell=False, index_list=None
    ):
        """
        Rotate atoms based on a vector and an angle, or two vectors. This function is completely adopted from ASE code
        (https://wiki.fysik.dtu.dk/ase/_modules/ase/atoms.html#Atoms.rotate)

        Args:

            rotate_cell:
            center:
            vector (list/numpy.ndarray/string):
                Vector to rotate the atoms around. Vectors can be given as
                strings: 'x', '-x', 'y', ... .

            angle (float/list) in radians = None:
                Angle that the atoms is rotated around the vecor 'v'. If an angle
                is not specified, the length of 'v' is used as the angle
                (default). The angle can also be a vector and then 'v' is rotated
                into 'a'.

            center = [0, 0, 0]:
                The center is kept fixed under the rotation. Use 'COM' to fix
                the center of mass, 'COP' to fix the center of positions or
                'COU' to fix the center of cell.

            rotate_cell = False:
                If true the cell is also rotated.

            index_list (list/numpy.ndarray):
                Indices of atoms to be rotated

        Examples:

        Rotate 90 degrees around the z-axis, so that the x-axis is
        rotated into the y-axis:

        >>> atoms = Atoms('H', [[-0.1, 1.01, -0.5]], cell=[[1, 0, 0], [0, 1, 0], [0, 0, 4]], pbc=[1, 1, 0])
        >>> a = (22./ 7.) / 2. # pi/2
        >>> atoms.rotate('z', a)
        >>> atoms.rotate((0, 0, 1), a)
        >>> atoms.rotate('-z', -a)
        >>> atoms.rotate((0, 0, a))
        >>> atoms.rotate('x', 'y')
        """

        norm = np.linalg.norm
        vector = string2vector(vector)
        if angle is None:
            angle = norm(vector)
        if isinstance(angle, (float, int)):
            vector /= norm(vector)
            c = cos(angle)
            s = sin(angle)
        else:
            v2 = string2vector(angle)
            vector /= norm(vector)
            v2 /= norm(v2)
            c = np.dot(vector, v2)
            vector = np.cross(vector, v2)
            s = norm(vector)
            # In case *v* and *a* are parallel, np.cross(v, v2) vanish
            # and can't be used as a rotation axis. However, in this
            # case any rotation axis perpendicular to v2 will do.
            eps = 1e-7
            if s < eps:
                vector = np.cross((0, 0, 1), v2)
                if norm(vector) < eps:
                    vector = np.cross((1, 0, 0), v2)
                if not (norm(vector) >= eps):
                    raise AssertionError()
            elif s > 0:
                vector /= s

        if isinstance(center, str):
            if center.lower() == "com":
                center = self.get_center_of_mass()
            elif center.lower() == "cop":
                center = np.mean(self.get_positions(), axis=0)
            elif center.lower() == "cou":
                center = self.cell.sum(axis=0) / 2
            else:
                raise ValueError("Cannot interpret center")
        else:
            center = np.array(center)

        if index_list is not None:
            if not (len(index_list) > 0):
                raise AssertionError()
            rotate_list = np.array(index_list)
        else:
            rotate_list = np.array(len(self) * [True])

        p = self.positions[rotate_list] - center
        self.positions[rotate_list] = (
            c * p
            - np.cross(p, s * vector)
            + np.outer(np.dot(p, vector), (1.0 - c) * vector)
            + center
        )
        if rotate_cell:
            rotcell = self.cell
            rotcell[:] = (
                c * rotcell
                - np.cross(rotcell, s * vector)
                + np.outer(np.dot(rotcell, vector), (1.0 - c) * vector)
            )
            self.set_cell(rotcell)

    def rotate_euler(self, center=(0, 0, 0), phi=0.0, theta=0.0, psi=0.0):
        """Rotate atoms via Euler angles.

        See e.g http://mathworld.wolfram.com/EulerAngles.html for explanation.

        Parameters:

        center :
            The point to rotate about. a sequence of length 3 with the
            coordinates, or 'COM' to select the center of mass, 'COP' to
            select center of positions or 'COU' to select center of cell.
        phi :
            The 1st rotation angle around the z axis (in radian)
        theta :
            Rotation around the x axis (in radian)
        psi :
            2nd rotation around the z axis (in radian)

        """
        if isinstance(center, str):
            if center.lower() == "com":
                center = self.get_center_of_mass()
            elif center.lower() == "cop":
                center = self.get_positions().mean(axis=0)
            elif center.lower() == "cou":
                center = self.cell.sum(axis=0) / 2
            else:
                raise ValueError("Cannot interpret center")
        else:
            center = np.array(center)

        # First move the molecule to the origin In contrast to MATLAB,
        # numpy broadcasts the smaller array to the larger row-wise,
        # so there is no need to play with the Kronecker product.
        if self._is_scaled:
            rcoords = self.get_scaled_positions(wrap=False) - center
        else:
            rcoords = self.positions - center

        # First Euler rotation about z in matrix form
        d = np.array(
            ((cos(phi), sin(phi), 0.0), (-sin(phi), cos(phi), 0.0), (0.0, 0.0, 1.0))
        )
        # Second Euler rotation about x:
        c = np.array(
            (
                (1.0, 0.0, 0.0),
                (0.0, cos(theta), sin(theta)),
                (0.0, -sin(theta), cos(theta)),
            )
        )
        # Third Euler rotation, 2nd rotation about z:
        b = np.array(
            ((cos(psi), sin(psi), 0.0), (-sin(psi), cos(psi), 0.0), (0.0, 0.0, 1.0))
        )
        # Total Euler rotation
        a = np.dot(b, np.dot(c, d))
        # Do the rotation
        rcoords = np.dot(a, np.transpose(rcoords))
        # Move back to the rotation point
        if self._is_scaled:
            self.set_scaled_positions(np.transpose(rcoords) + center)
        else:
            self.positions = np.transpose(rcoords) + center

    def set_calculator(self, calc=None):
        """Attach calculator object."""
        pass
     
    def get_calculator(self):
        """Get currently attached calculator object."""
        return None
        
    def translate(self, displacement):
        """
        Translate atomic positions.

        The displacement argument can be a float, an xyz vector, or an
        nx3 array (where n is the number of atoms).

        Args:
            displacement:

        Returns:

        """
        self.positions += np.array(displacement)

    def wrap(self, center=(0.5, 0.5, 0.5), pbc=None, eps=1e-7):
        """Wrap positions to unit cell.

        Parameters:

        center: three float
            The positons in fractional coordinates that the new positions
            will be nearest possible to.
        pbc: one or 3 bool
            For each axis in the unit cell decides whether the positions
            will be moved along this axis.  By default, the boundary
            conditions of the Atoms object will be used.
        eps: float
            Small number to prevent slightly negative coordinates from beeing
            wrapped.

        See also the :func:`ase.utils.geometry.wrap_positions` function.
        Example:

        >>> a = Atoms('H',
        ...           [[-0.1, 1.01, -0.5]],
        ...           cell=[[1, 0, 0], [0, 1, 0], [0, 0, 4]],
        ...           pbc=[1, 1, 0])
        >>> a.wrap()
        >>> a.positions
        array([[ 0.9 ,  0.01, -0.5 ]])
        """

        from ase.utils.geometry import wrap_positions

        if pbc is None:
            pbc = self.pbc
        self.positions = wrap_positions(
            positions=self.positions,
            cell=self.cell,
            pbc=pbc,
            center=center,
            eps=eps
        )

    def write(self, filename, format=None, **kwargs):
        """
        Write atoms object to a file.

        see ase.io.write for formats.
        kwargs are passed to ase.io.write.

        Args:
            filename:
            format:
            **kwargs:

        Returns:

        """
        from ase.io import write

        atoms = self.copy()
        atoms.arrays["positions"] = atoms.positions
        write(filename, atoms, format, **kwargs)


class _CrystalStructure(Atoms):
    """
    only for historical reasons

    Args:
        element:
        BravaisLattice:
        BravaisBasis:
        LatticeConstants:
        Dimension:
        relCoords:
        PSE:
        **kwargs:
    """

    def __init__(
        self,
        element="Fe",
        bravais_lattice="cubic",
        bravais_basis="primitive",
        lattice_constants=None,  # depending on symmetry length and angles
        dimension=3,
        rel_coords=True,
        pse=None,
        **kwargs
    ):

        # print "basis0"
        # allow also for scalar input for LatticeConstants (for a cubic system)
        if lattice_constants is None:
            lattice_constants = [1.0]
        try:
            test = lattice_constants[0]
        except (TypeError, IndexError):
            lattice_constants = [lattice_constants]
        self.bravais_lattice = bravais_lattice
        self.bravais_basis = bravais_basis
        self.lattice_constants = lattice_constants
        self.dimension = dimension
        self.relCoords = rel_coords
        self.element = element

        self.__updateCrystal__(pse)

        self.crystalParamsDict = {
            "BravaisLattice": self.bravais_lattice,
            "BravaisBasis": self.bravais_basis,
            "LatticeConstants": self.lattice_constants,
        }

        self.crystal_lattice_dict = {
            3: {
                "cubic": ["fcc", "bcc", "primitive"],
                "hexagonal": ["primitive", "hcp"],
                "monoclinic": ["primitive", "base-centered"],
                "triclinic": ["primitive"],
                "orthorombic": [
                    "primitive",
                    "body-centered",
                    "base-centered",
                    "face-centered",
                ],
                "tetragonal": ["primitive", "body-centered"],
                "rhombohedral": ["primitive"],
            },
            2: {
                "oblique": ["primitive"],
                "rectangular": ["primitive", "centered"],
                "hexagonal": ["primitive"],
                "square": ["primitive"],
            },
            1: {"line": ["primitive"]},
        }

        # init structure for lattice parameters alat, blat, clat, alpha, beta, gamma
        self.crystalLatticeParams = {
            3: {
                "cubic": [1.0],
                "hexagonal": [1.0, 2.0],
                "monoclinic": [1.0, 1.0, 1.0, 90.0],
                "triclinic": [1.0, 2.0, 3.0, 90.0, 90.0, 90.0],
                "orthorombic": [1.0, 1.0, 1.0],
                "tetragonal": [1.0, 2.0],
                "rhombohedral": [1.0, 90.0, 90.0, 90.0],
            },
            2: {
                "oblique": [1.0, 1.0, 90.0],
                "rectangular": [1.0, 1.0],
                "hexagonal": [1.0],
                "square": [1.0],
            },
            1: {"line": [1.0]},
        }

        # print "basis"
        super(_CrystalStructure, self).__init__(
            elements=self.ElementList,
            scaled_positions=self.coordinates,
            cell=self.amat,  # tag = "Crystal",
            pbc=[True, True, True][0 : self.dimension],
        )

    # ## private member functions
    def __updateCrystal__(self, pse=None):
        """

        Args:
            pse:

        Returns:

        """
        self.__updateAmat__()
        self.__updateCoordinates__()
        self.__updateElementList__(pse)

    def __updateAmat__(self):  # TODO: avoid multi-call of this function
        """

        Returns:

        """
        # print "lat constants (__updateAmat__):", self.LatticeConstants
        a_lat = self.lattice_constants[0]

        if self.dimension == 3:
            alpha = None
            beta = None
            gamma = None
            b_lat, c_lat = None, None
            if self.bravais_lattice == "cubic":
                b_lat = c_lat = a_lat
                alpha = beta = gamma = 90 / 180.0 * np.pi  # 90 degrees
            elif self.bravais_lattice == "tetragonal":
                b_lat = a_lat
                c_lat = self.lattice_constants[1]
                alpha = beta = gamma = 0.5 * np.pi  # 90 degrees
            elif self.bravais_lattice == "triclinic":
                b_lat = self.lattice_constants[1]
                c_lat = self.lattice_constants[2]
                alpha = self.lattice_constants[3] / 180.0 * np.pi
                beta = self.lattice_constants[4] / 180.0 * np.pi
                gamma = self.lattice_constants[5] / 180.0 * np.pi
            elif self.bravais_lattice == "hexagonal":
                b_lat = a_lat
                c_lat = self.lattice_constants[1]
                alpha = 60.0 / 180.0 * np.pi  # 60 degrees
                beta = gamma = 0.5 * np.pi  # 90 degrees
            elif self.bravais_lattice == "orthorombic":
                b_lat = self.lattice_constants[1]
                c_lat = self.lattice_constants[2]
                alpha = beta = gamma = 0.5 * np.pi  # 90 degrees
            elif self.bravais_lattice == "rhombohedral":
                b_lat = a_lat
                c_lat = a_lat
                alpha = self.lattice_constants[1] / 180.0 * np.pi
                beta = self.lattice_constants[2] / 180.0 * np.pi
                gamma = self.lattice_constants[3] / 180.0 * np.pi
            elif self.bravais_lattice == "monoclinic":
                b_lat = self.lattice_constants[1]
                c_lat = self.lattice_constants[2]
                alpha = 0.5 * np.pi
                beta = self.lattice_constants[3] / 180.0 * np.pi
                gamma = 0.5 * np.pi

            b1 = np.cos(alpha)
            b2 = np.sin(alpha)
            c1 = np.cos(beta)
            c2 = (np.cos(gamma) - np.cos(beta) * np.cos(alpha)) / np.sin(alpha)
            self.amat = np.array(
                [
                    [a_lat, 0.0, 0.0],
                    [b_lat * b1, b_lat * b2, 0.0],
                    [c_lat * c1, c_lat * c2, c_lat * np.sqrt(1 - c2 * c2 - c1 * c1)],
                ]
            )
        elif self.dimension == 2:  # TODO not finished yet
            self.amat = a_lat * np.array([[1.0, 0.0], [0.0, 1.0]])
            if self.bravais_lattice == "rectangular":
                b_lat = self.lattice_constants[1]
                self.amat = np.array([[a_lat, 0.0], [0.0, b_lat]])
        elif self.dimension == 1:
            self.amat = a_lat * np.array([[1.0]])
        else:
            raise ValueError("Bravais lattice not defined!")

    def __updateElementList__(self, pse=None):
        """

        Args:
            pse:

        Returns:

        """
        self.ElementList = len(self.coordinates) * [self.element]

    def __updateCoordinates__(self):
        """

        Returns:

        """
        # if relative coordinates
        basis = None
        if self.dimension == 3:
            if self.bravais_basis == "fcc" or self.bravais_basis == "face-centered":
                basis = np.array(
                    [[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]]
                )
            elif self.bravais_basis == "body-centered" or self.bravais_basis == "bcc":
                basis = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
            elif self.bravais_basis == "base-centered":
                basis = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0]])
            elif self.bravais_basis == "hcp":
                # basis = r([[0.0,-1./np.sqrt(3.),np.sqrt(8./3.)]])
                # a = self.LatticeConstants[0]
                # c = self.LatticeConstants[1]
                basis = np.array([[0.0, 0.0, 0.0], [1.0 / 3.0, 1.0 / 3.0, 1.0 / 2.0]])
                # basis = np.dot(basis,np.linalg.inv(self.amat))
            elif self.bravais_basis == "primitive":
                basis = np.array([[0.0, 0.0, 0.0]])
            else:
                exit()
        elif self.dimension == 2:
            if self.bravais_basis == "primitive":
                basis = np.array([[0.0, 0.0]])
            elif self.bravais_basis == "centered":
                basis = np.array([[0.0, 0.0], [0.5, 0.5]])
            else:
                exit()
        elif self.dimension == 1:
            if self.bravais_basis == "primitive":
                basis = np.array([[0.0]])
            else:
                exit()
        self.coordinates = basis

    # ########################### get commmands ########################
    def get_lattice_types(self):
        """

        Returns:

        """
        self.crystal_lattice_dict[self.dimension].keys().sort()
        return self.crystal_lattice_dict[self.dimension].keys()

    def get_dimension_of_lattice_parameters(self):
        """

        Returns:

        """
        # print "getDimensionOfLatticeParameters"
        counter = 0
        for k in self.get_needed_lattice_parameters():
            if k:
                counter += 1
        return counter

    def get_needed_lattice_parameters(self):
        """

        Returns:

        """
        # print "call: getNeededLatticeParams"
        needed_params = [True, False, False, False, False, False]
        if self.dimension == 3:
            if self.bravais_lattice == "cubic":
                needed_params = [
                    True,
                    False,
                    False,
                    False,
                    False,
                    False,
                ]  # stands for alat, blat, clat, alpha, beta, gamma
            elif self.bravais_lattice == "triclinic":
                needed_params = [True, True, True, True, True, True]
            elif self.bravais_lattice == "monoclinic":
                needed_params = [True, True, True, True, False, False]
            elif self.bravais_lattice == "orthorombic":
                needed_params = [True, True, True, False, False, False]
            elif self.bravais_lattice == "tetragonal":
                needed_params = [True, False, True, False, False, False]
            elif self.bravais_lattice == "rhombohedral":
                needed_params = [True, False, False, True, True, True]
            elif self.bravais_lattice == "hexagonal":
                needed_params = [True, False, True, False, False, False]
        elif self.dimension == 2:
            if self.bravais_lattice == "oblique":
                needed_params = [True, True, False, True, False, False]
            elif self.bravais_lattice == "rectangular":
                needed_params = [True, True, False, False, False, False]
            elif self.bravais_lattice == "hexagonal":
                needed_params = [True, False, False, False, False, False]
            elif self.bravais_lattice == "square":
                needed_params = [True, False, False, False, False, False]
            else:  # TODO: need to be improved
                needed_params = [True, False, False, False, False, False]
        elif self.dimension == 1:
            if self.bravais_lattice == "line":
                needed_params = [True, False, False, False, False, False]
            else:  # TODO: improval needed
                needed_params = [True, False, False, False, False, False]
        else:
            raise ValueError("inconsistency in lattice structures")

        return needed_params

    def get_basis_types(self):
        """

        Returns:

        """
        self.crystal_lattice_dict[self.dimension].get(self.bravais_lattice).sort()
        return self.crystal_lattice_dict[self.dimension].get(self.bravais_lattice)

    def get_initial_lattice_constants(self):
        """

        Returns:

        """
        self.crystalLatticeParams[self.dimension].get(self.bravais_lattice).sort()
        return (
            self.crystalLatticeParams[self.dimension].get(self.bravais_lattice).sort()
        )

    # def getDimension(self):
    #     return self.dimension

    # def getCoordinates(self):
    #     return self.coordinates

    # def getCell(self):
    #     return self.amat

    def get_atom_structure(self, rel=True):
        """

        Args:
            rel:

        Returns:

        """
        #        print self.relCoords, self.amat
        return Atoms(
            elementList=self.ElementList,
            coordinates=self.coordinates,
            amat=self.amat,
            tag="Crystal",
            rel=rel,  # self.relCoords, #rel, # true or false # coordinates are given in relative lattice units
            pbc=[True, True, True][0 : self.dimension],
            Crystal=self.crystalParamsDict,
        )

    # #################### set commands #########################
    def set_lattice_constants(self, lattice_constants=None):
        """

        Args:
            lattice_constants:

        Returns:

        """
        if lattice_constants is None:
            lattice_constants = [1.0]
        for k in lattice_constants:
            if k <= 0:
                raise ValueError("negative lattice parameter(s)")
        self.lattice_constants = lattice_constants
        self.__updateCrystal__()

    def set_element(self, element="Fe"):
        """

        Args:
            element:

        Returns:

        """
        self.element = element
        self.__updateCrystal__()

    def set_dimension(self, dim=3):
        """

        Args:
            dim:

        Returns:

        """
        self.dimension = dim
        length = self.get_dimension_of_lattice_parameters()
        if dim == 3:  # # initial 3d structure
            self.lattice_constants = length * [1.0]
            self.bravais_lattice = "cubic"
            self.bravais_basis = "primitive"
        elif dim == 2:  # # initial 2d structure
            self.lattice_constants = length * [1.0]
            self.bravais_lattice = "square"
            self.bravais_basis = "primitive"
        elif dim == 1:  # # initial 1d structure
            self.lattice_constants = length * [1.0]
            self.bravais_lattice = "line"
            self.bravais_basis = "primitive"
        self.__updateCrystal__()

    def set_lattice_type(self, name_lattice="cubic"):
        """

        Args:
            name_lattice:

        Returns:

        """
        # catch input error
        # print "lattice type =", name_lattice
        if name_lattice not in self.get_lattice_types():
            raise ValueError("is not item of ")
        else:
            self.bravais_lattice = name_lattice
            self.set_lattice_constants(
                self.get_dimension_of_lattice_parameters() * [1.0]
            )
            self.set_basis_type(
                name_basis=self.crystal_lattice_dict[self.dimension].get(name_lattice)[
                    0
                ]
            )  # initial basis type

        self.__updateCrystal__()

    def set_basis_type(self, name_basis="primitive"):
        """

        Args:
            name_basis:

        Returns:

        """
        if name_basis not in self.get_basis_types():
            raise ValueError("is not item of")
        else:
            self.bravais_basis = name_basis
        self.__updateCrystal__()

    def atoms(self):
        """

        Returns:

        """
        return Atoms(
            elements=self.ElementList,
            scaled_positions=self.coordinates,
            cell=self.amat,
            pbc=[True, True, True][0 : self.dimension],
        )


class Neighbors:
    """
    Class for storage of the neighbor information for a given atom based on the KDtree algorithm
    """

    def __init__(self):
        self._distances = None
        self._vecs = None
        self._indices = None
        self._shells = None

    @property
    def distances(self):
        return self._distances

    @distances.setter
    def distances(self, new_distances):
        if isinstance(new_distances, list) or isinstance(new_distances, np.ndarray):
            self._distances = np.array(new_distances)
        else:
            raise TypeError("Only lists and np.arrays are supported.")

    @property
    def vecs(self):
        return self._vecs

    @vecs.setter
    def vecs(self, new_vecs):
        if isinstance(new_vecs, list) or isinstance(new_vecs, np.ndarray):
            self._vecs = np.array(new_vecs)
        else:
            raise TypeError("Only lists and np.arrays are supported.")

    @property
    def indices(self):
        return self._indices

    @indices.setter
    def indices(self, new_indices):
        if isinstance(new_indices, list) or isinstance(new_indices, np.ndarray):
            self._indices = np.array(new_indices)
        else:
            raise TypeError("Only lists and np.arrays are supported.")

    @property
    def shells(self):
        return self._shells

    @shells.setter
    def shells(self, new_shells):
        if isinstance(new_shells, list) or isinstance(new_shells, np.array):
            self._shells = np.array(new_shells)
        else:
            raise TypeError("Only lists and np.arrays are supported.")


class CrystalStructure(object):
    def __new__(cls, *args, **kwargs):
        basis = _CrystalStructure(*args, **kwargs).atoms()
        return basis


def ase_to_pyiron(ase_obj):
    """

    Args:
        ase_obj:

    Returns:

    """
    try:
        import ase
    except ImportError:
        raise ValueError("ASE package not yet installed")
    element_list = ase_obj.get_chemical_symbols()
    cell = ase_obj.cell
    positions = ase_obj.get_positions()
    pbc = ase_obj.get_pbc()
    spins = ase_obj.get_initial_magnetic_moments()
    if all(spins == np.array(None)) or sum(np.abs(spins)) == 0.0:
        pyiron_atoms = Atoms(
            elements=element_list, positions=positions, pbc=pbc, cell=cell
        )
    else:
        if any(spins == np.array(None)):
            spins[spins == np.array(None)] = 0.0
        pyiron_atoms = Atoms(
            elements=element_list,
            positions=positions,
            pbc=pbc,
            cell=cell,
            magmoms=spins,
        )
    if hasattr(ase_obj, "constraints") and len(ase_obj.constraints) != 0:
        for constraint in ase_obj.constraints:
            constraint_dict = constraint.todict()
            if constraint_dict["name"] == "FixAtoms":
                if "selective_dynamics" not in pyiron_atoms._tag_list.keys():
                    pyiron_atoms.add_tag(selective_dynamics=[True, True, True])
                pyiron_atoms.selective_dynamics[
                    constraint_dict["kwargs"]["indices"]
                ] = [False, False, False]
            elif constraint_dict["name"] == "FixScaled":
                if "selective_dynamics" not in pyiron_atoms._tag_list.keys():
                    pyiron_atoms.add_tag(selective_dynamics=[True, True, True])
                pyiron_atoms.selective_dynamics[
                    constraint_dict["kwargs"]["a"]
                ] = constraint_dict["kwargs"]["mask"]
            else:
                warnings.warn("Unsupported ASE constraint: " + constraint_dict["name"])
    return pyiron_atoms


def pyiron_to_ase(pyiron_obj):
    try:
        from pyiron.atomistics.structure.pyironase import ASEAtoms
    except ImportError:
        raise ValueError("ASE package not yet installed")
    element_list = pyiron_obj.get_parent_symbols()
    cell = pyiron_obj.cell
    positions = pyiron_obj.positions
    pbc = pyiron_obj.get_pbc()
    spins = pyiron_obj.get_initial_magnetic_moments()
    if all(spins == np.array(None)) or sum(np.abs(spins)) == 0.0:
        atoms = ASEAtoms(symbols=element_list, positions=positions, pbc=pbc, cell=cell)
    else:
        if any(spins == np.array(None)):
            spins[spins == np.array(None)] = 0.0
        atoms = ASEAtoms(
            symbols=element_list, positions=positions, pbc=pbc, cell=cell, magmoms=spins
        )
    return atoms


def _check_if_simple_atoms(atoms):
    """
    Raise a warning if the ASE atoms object includes properties which can not be converted to pymatgen atoms.

    Args:
        atoms: ASE atoms object
    """
    dict_keys = [
        k
        for k in atoms.__dict__.keys()
        if k
        not in ["_celldisp", "arrays", "_cell", "_pbc", "_constraints", "info", "_calc"]
    ]
    array_keys = [
        k for k in atoms.__dict__["arrays"].keys() if k not in ["numbers", "positions"]
    ]
    if not len(dict_keys) == 0:
        warnings.warn("Found unknown keys: " + str(dict_keys))
    if not np.all(atoms.__dict__["_celldisp"] == np.array([[0.0], [0.0], [0.0]])):
        warnings.warn("Found cell displacement: " + str(atoms.__dict__["_celldisp"]))
    if not atoms.__dict__["_calc"] is None:
        warnings.warn("Found calculator: " + str(atoms.__dict__["_calc"]))
    if not atoms.__dict__["_constraints"] == []:
        warnings.warn("Found constraint: " + str(atoms.__dict__["_constraints"]))
    if not np.all(atoms.__dict__["_pbc"]):
        warnings.warn("Cell is not periodic: " + str(atoms.__dict__["_pbc"]))
    if not len(array_keys) == 0:
        warnings.warn("Found unknown flags: " + str(array_keys))
    if not atoms.__dict__["info"] == dict():
        warnings.warn("Info is not empty: " + str(atoms.__dict__["info"]))


def pymatgen_to_pyiron(pymatgen_obj):
    """
    Convert pymatgen atoms object to pyiron atoms object (pymatgen->ASE->pyiron)

    Args:
        pymatgen_obj: pymatgen atoms object

    Returns:
        pyiron atoms object
    """
    try:
        from pymatgen.io.ase import AseAtomsAdaptor
    except ImportError:
        raise ValueError("PyMatGen package not yet installed")
    return ase_to_pyiron(AseAtomsAdaptor().get_atoms(structure=pymatgen_obj))


def pyiron_to_pymatgen(pyiron_obj):
    """
    Convert pyiron atoms object to pymatgen atoms object

    Args:
        pyiron_obj: pyiron atoms object

    Returns:
        pymatgen atoms object
    """
    try:
        from pymatgen.io.ase import AseAtomsAdaptor
    except ImportError:
        raise ValueError("PyMatGen package not yet installed")
    ase_atoms = pyiron_to_ase(pyiron_obj)
    _check_if_simple_atoms(atoms=ase_atoms)
    return AseAtomsAdaptor().get_structure(atoms=ase_atoms, cls=None)


def ovito_to_pyiron(ovito_obj):
    """

    Args:
        ovito_obj:

    Returns:

    """
    try:
        from ovito.data import ase_to_pyiron

        return ase_to_pyiron(ovito_obj.to_ase_atoms())
    except ImportError:
        raise ValueError("ovito package not yet installed")


def pyiron_to_ovito(atoms):
    """

    Args:
        atoms:

    Returns:

    """
    try:
        from ovito.data import DataCollection

        return DataCollection.create_from_ase_atoms(atoms)
    except ImportError:
        raise ValueError("ovito package not yet installed")


def string2symbols(s):
    """
    Convert string to list of chemical symbols.

    Args:
        s:

    Returns:

    """
    i = None
    n = len(s)
    if n == 0:
        return []
    c = s[0]
    if c.isdigit():
        i = 1
        while i < n and s[i].isdigit():
            i += 1
        return int(s[:i]) * string2symbols(s[i:])
    if c == "(":
        p = 0
        for i, c in enumerate(s):
            if c == "(":
                p += 1
            elif c == ")":
                p -= 1
                if p == 0:
                    break
        j = i + 1
        while j < n and s[j].isdigit():
            j += 1
        if j > i + 1:
            m = int(s[i + 1 : j])
        else:
            m = 1
        return m * string2symbols(s[1:i]) + string2symbols(s[j:])

    if c.isupper():
        i = 1
        if 1 < n and s[1].islower():
            i += 1
        j = i
        while j < n and s[j].isdigit():
            j += 1
        if j > i:
            m = int(s[i:j])
        else:
            m = 1
        return m * [s[:i]] + string2symbols(s[j:])
    else:
        raise ValueError


def symbols2numbers(symbols):
    """

    Args:
        symbols (list, str):

    Returns:

    """
    pse = PeriodicTable()
    df = pse.dataframe.T
    if isinstance(symbols, str):
        symbols = string2symbols(symbols)
    numbers = list()
    for sym in symbols:
        if isinstance(sym, string_types):
            numbers.append(df[sym]["AtomicNumber"])
        else:
            numbers.append(sym)
    return numbers


def string2vector(v):
    """

    Args:
        v:

    Returns:

    """
    if isinstance(v, str):
        if v[0] == "-":
            return -string2vector(v[1:])
        w = np.zeros(3)
        w["xyz".index(v)] = 1.0
        return w
    return np.array(v, float)


def default(data, dflt):
    """
    Helper function for setting default values.

    Args:
        data:
        dflt:

    Returns:

    """
    if data is None:
        return None
    elif isinstance(data, (list, tuple)):
        newdata = []
        allnone = True
        for x in data:
            if x is None:
                newdata.append(dflt)
            else:
                newdata.append(x)
                allnone = False
        if allnone:
            return None
        return newdata
    else:
        return data
