# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
from collections import OrderedDict
import numpy as np
from pyiron_base import GenericParameters
import decimal as dec

try:
    from ase.calculators.lammps import Prism

except ImportError:
    try:
        from ase.calculators.lammpsrun import Prism
    except ImportError:
        from ase.calculators.lammpsrun import prism as Prism

__author__ = "Joerg Neugebauer, Sudarsan Surendralal, Yury Lysogorskiy, Jan Janssen, Markus Tautschnig"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Sudarsan Surendralal"
__email__ = "surendralal@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


class UnfoldingPrism(Prism):
    """
    Create a lammps-style triclinic prism object from a cell

    The main purpose of the prism-object is to create suitable
    string representations of prism limits and atom positions
    within the prism.
    When creating the object, the digits parameter (default set to 10)
    specify the precision to use.
    lammps is picky about stuff being within semi-open intervals,
    e.g. for atom positions (when using create_atom in the in-file),
    x must be within [xlo, xhi).

    Args:
        cell:
        pbc:
        digits:
    """

    def __init__(self, cell, pbc=(True, True, True), digits=10):
        # Temporary fix. Since the arguments for the constructor have changed, try to see if it is compatible with
        # the latest ase. If not, revert to the old __init__ parameters.
        try:
            super(UnfoldingPrism, self).__init__(
                cell, pbc=pbc, tolerance=float("1e-{}".format(digits))
            )
        except TypeError:
            super(UnfoldingPrism, self).__init__(cell, pbc=pbc, digits=digits)
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
        self.car_prec = dec.Decimal("10.0") ** int(
            np.floor(np.log10(max((xhi, yhi, zhi)))) - digits
        )
        self.dir_prec = dec.Decimal("10.0") ** (-digits)
        self.acc = float(self.car_prec)
        self.eps = np.finfo(xhi).eps

        # For rotating positions from ase to lammps
        apre = np.array(((xhi, 0, 0), (xyp, yhi, 0), (xzp, yzp, zhi)))
        # np.linalg.inv(cell) ?= np.array([np.cross(b, c), np.cross(c, a), np.cross(a, b)]).T / np.linalg.det(cell)
        self.R = np.dot(np.linalg.inv(cell), apre)

        def fold(vec, pvec, i):
            p = pvec[i]
            x = vec[i] + 0.5 * p
            n = (np.mod(x, p) - x) / p
            return [float(self.f2qdec(vec_a)) for vec_a in (vec + n * pvec)], n

        apre[1, :], n1 = fold(apre[1, :], apre[0, :], 0)
        if np.abs(apre[1, 0] / apre[0, 0]) > 0.5:
            apre[1, 0] -= np.sign(n1) * apre[0, 0]
            n1 -= np.sign(n1)

        apre[2, :], n2 = fold(apre[2, :], apre[1, :], 1)
        if np.abs(apre[2, 1] / apre[1, 1]) > 0.5:
            apre[2, 1] -= np.sign(n2) * apre[1, 1]
            n2 -= np.sign(n2)

        apre[2, :], n3 = fold(apre[2, :], apre[0, :], 0)
        if np.abs(apre[2, 0] / apre[0, 0]) > 0.5:
            apre[2, 0] -= np.sign(n3) * apre[0, 0]
            n3 -= np.sign(n3)
        self.ns = [n1, n2, n3]

        d_a = apre[0, 0] / 2 - apre[1, 0]
        if np.abs(d_a) < self.acc:
            if d_a < 0:
                print("debug: apply shift")
                apre[1, 0] += 2 * d_a
                apre[2, 0] += 2 * d_a

        self.A = apre

        if self.is_skewed() and (not (pbc[0] and pbc[1] and pbc[2])):
            raise RuntimeError(
                "Skewed lammps cells MUST have " "PBC == True in all directions!"
            )

    def unfold_cell(self, cell):
        """
        Unfold LAMMPS cell to original

        Let C be the pyiron cell and A be the Lammps cell, then define (in init) the rotation matrix between them as
            R := C^inv.A
        And recall that rotation matrices have the property
            R^T == R^inv
        Then left multiply the definition of R by C, and right multiply by R.T to get
            C.R.R^T = C.C^inv.A.R^T
        Then
            C = A.R^T

        After that, account for the folding process.

        Args:
            cell: LAMMPS cell,

        Returns:
            unfolded cell
        """
        # Rotation
        ucell = np.dot(cell, self.R.T)
        # Folding
        a = ucell[0]
        bp = ucell[1]
        cpp = ucell[2]
        (n1, n2, n3) = self.ns
        b = bp - n1 * a
        c = cpp - n2 * bp - n3 * a
        return np.array([a, b, c])

    def pos_to_lammps(self, position):
        """
        Rotate an ase-cell position to the lammps cell orientation

        Args:
            position:

        Returns:
            tuple of float.
        """
        return tuple([x for x in np.dot(position, self.R)])

    def f2qdec(self, f):
        return dec.Decimal(repr(f)).quantize(self.car_prec, dec.ROUND_DOWN)

    def f2s(self, f):
        return str(dec.Decimal(repr(f)).quantize(self.car_prec, dec.ROUND_HALF_EVEN))

    def get_lammps_prism_str(self):
        """Return a tuple of strings"""
        p = self.get_lammps_prism()
        return tuple([self.f2s(x) for x in p])


class LammpsStructure(GenericParameters):
    """

    Args:
        input_file_name:
    """

    def __init__(self, input_file_name=None, bond_dict=None):
        super(LammpsStructure, self).__init__(
            input_file_name=input_file_name,
            table_name="structure_inp",
            comment_char="#",
            val_only=True,
        )
        self._structure = None
        self._potential = None
        self._el_eam_lst = []
        self.atom_type = None
        self.cutoff_radius = None
        self.digits = 10
        self._bond_dict = bond_dict

    @property
    def potential(self):
        return self._potential

    @potential.setter
    def potential(self, val):
        self._potential = val

    @property
    def structure(self):
        """

        Returns:

        """
        return self._structure

    @structure.setter
    def structure(self, structure):
        """

        Args:
            structure:

        Returns:

        """
        self._structure = structure
        if self.atom_type == "full":
            input_str = self.structure_full()
        elif self.atom_type == "bond":
            input_str = self.structure_bond()
        elif self.atom_type == "charge":
            input_str = self.structure_charge()
        else:  # self.atom_type == 'atomic'
            input_str = self.structure_atomic()
        self.load_string(input_str)

    @property
    def el_eam_lst(self):
        """

        Returns:

        """
        return self._el_eam_lst

    @el_eam_lst.setter
    def el_eam_lst(self, el_eam_lst):
        """

        Args:
            el_eam_lst:

        Returns:

        """
        self._el_eam_lst = el_eam_lst

    def load_default(self):
        """

        Returns:

        """
        input_str = ""
        self.load_string(input_str)

    def simulation_cell(self):
        """

        Returns:

        """

        self.prism = UnfoldingPrism(self._structure.cell, digits=15)
        xhi, yhi, zhi, xy, xz, yz = self.prism.get_lammps_prism_str()
        # Please, be carefull and not round xhi, yhi,..., otherwise you will get too skew cell from LAMMPS.
        # These values are already checked in UnfoldingPrism to fullfill LAMMPS skewness criteria
        simulation_cell = (
            "0. {} xlo xhi\n".format(xhi)
            + "0. {} ylo yhi\n".format(yhi)
            + "0. {} zlo zhi\n".format(zhi)
        )

        if self.prism.is_skewed():
            simulation_cell += "{0} {1} {2} xy xz yz\n".format(xy, xz, yz)

        return simulation_cell

    def structure_bond(self):
        """

        Returns:

        """
        # analyze structure to get molecule_ids, bonds, angles etc
        coords = self.rotate_positions(self._structure)

        elements = self._structure.get_chemical_symbols()
        el_list = self._structure.get_species_symbols()
        el_dict = OrderedDict()
        for object_id, el in enumerate(el_list):
            el_dict[el] = object_id

        n_s = len(el_list)
        bond_type = np.ones([n_s, n_s], dtype=np.int)
        count = 0
        for i in range(n_s):
            for j in range(i, n_s):
                count += 1
                bond_type[i, j] = count
                bond_type[j, i] = count

        if self.structure.bonds is None:
            if self.cutoff_radius is None:
                bonds_lst = self.structure.get_bonds(max_shells=1)
            else:
                bonds_lst = self.structure.get_bonds(radius=self.cutoff_radius)
            bonds = []

            for ia, i_bonds in enumerate(bonds_lst):
                el_i = el_dict[elements[ia]]
                for el_j, b_lst in i_bonds.items():
                    b_type = bond_type[el_i][el_dict[el_j]]
                    for i_shell, ib_shell_lst in enumerate(b_lst):
                        for ib in np.unique(ib_shell_lst):
                            if ia < ib:  # avoid double counting of bonds
                                bonds.append([ia + 1, ib + 1, b_type])

            self.structure.bonds = np.array(bonds)
        bonds = self.structure.bonds

        atomtypes = (
            " Start File for LAMMPS \n"
            + "{0:d} atoms".format(len(self._structure))
            + " \n"
            + "{0:d} bonds".format(len(bonds))
            + " \n"
            + "{0} atom types".format(self._structure.get_number_of_species())
            + " \n"
            + "{0} bond types".format(np.max(bond_type))
            + " \n"
        )

        cell_dimesions = self.simulation_cell()

        masses = "Masses \n\n"
        el_obj_list = self._structure.get_species_objects()
        for object_id, el in enumerate(el_obj_list):
            masses += "{0:3d} {1:f}".format(object_id + 1, el.AtomicMass) + "\n"

        atoms = "Atoms \n\n"

        # atom_style bond
        # format: atom-ID, molecule-ID, atom_type, x, y, z
        format_str = "{0:d} {1:d} {2:d} {3:f} {4:f} {5:f} "
        if self._structure.dimension == 3:
            for id_atom, (x, y, z) in enumerate(coords):
                id_mol, id_species = 1, el_dict[elements[id_atom]]
                atoms += (
                    format_str.format(id_atom + 1, id_mol, id_species + 1, x, y, z)
                    + "\n"
                )
        elif self._structure.dimension == 2:
            for id_atom, (x, y) in enumerate(coords):
                id_mol, id_species = 1, el_dict[elements[id_atom]]
                atoms += (
                    format_str.format(id_atom + 1, id_mol, id_species + 1, x, y, 0.0)
                    + "\n"
                )
        else:
            raise ValueError("dimension 1 not yet implemented")

        bonds_str = "Bonds \n\n"
        for i_bond, (i_a, i_b, b_type) in enumerate(bonds):
            bonds_str += (
                "{0:d} {1:d} {2:d} {3:d}".format(i_bond + 1, b_type, i_a, i_b) + "\n"
            )

        return (
            atomtypes
            + "\n"
            + cell_dimesions
            + "\n"
            + masses
            + "\n"
            + atoms
            + "\n"
            + bonds_str
            + "\n"
        )

    def structure_full(self):
        """
        Write routine to create atom structure static file for atom_type='full' that can be loaded by LAMMPS

        Returns:

        """
        coords = self.rotate_positions(self._structure)

        # extract electric charges from potential file
        q_dict = {}
        species_translate_list = []
        for species in self.structure.species:
            species_name = species.Abbreviation
            q_dict[species_name] = self.potential.get_charge(species_name)
            species_translate_list.append(self.potential.get_element_id(species_name))
        sorted_species_list = np.array(self._potential.get_element_lst())
        molecule_lst, bonds_lst, angles_lst = [], [], []
        bond_type_lst, angle_type_lst = [], []
        # Using a cutoff distance to draw the bonds instead of the number of neighbors
        cutoff_list = list()
        for val in self._bond_dict.values():
            cutoff_list.append(np.max(val["cutoff_list"]))
        max_cutoff = np.max(cutoff_list)
        # Calculate neighbors only once
        neighbors = self._structure.get_neighbors_by_distance(cutoff_radius=max_cutoff)
        id_mol = 0
        indices = self._structure.indices
        for id_el, id_species in enumerate(indices):
            id_mol += 1
            molecule_lst.append([id_el, id_mol, species_translate_list[id_species]])
        # Draw bonds between atoms is defined in self._bond_dict
        # Go through all elements for which bonds are defined
        for element, val in self._bond_dict.items():
            el_1_list = self._structure.select_index(element)
            if el_1_list is not None:
                if len(el_1_list) > 0:
                    for i, v in enumerate(val["element_list"]):
                        el_2_list = self._structure.select_index(v)
                        cutoff_dist = val["cutoff_list"][i]
                        for j, ind in enumerate(np.array(neighbors.indices[el_1_list])):
                            # Only chose those indices within the cutoff distance and which belong
                            # to the species defined in the element_list
                            # i is the index of each bond type, and j is the element index
                            id_el = el_1_list[j]
                            bool_1 = np.array(neighbors.distances[id_el]) <= cutoff_dist
                            act_ind = ind[bool_1]
                            bool_2 = np.in1d(act_ind, el_2_list)
                            final_ind = act_ind[bool_2]
                            # Get the bond and angle type
                            bond_type = val["bond_type_list"][i]
                            angle_type = val["angle_type_list"][i]
                            # Draw only maximum allowed bonds
                            final_ind = final_ind[:val["max_bond_list"][i]]
                            for fi in final_ind:
                                bonds_lst.append([id_el + 1, fi + 1])
                                bond_type_lst.append(bond_type)
                            # Draw angles if at least 2 bonds are present and if an angle type is defined for this
                            # particular set of bonds
                            if len(final_ind) >= 2 and val["angle_type_list"][i] is not None:
                                angles_lst.append([final_ind[0] + 1, id_el + 1, final_ind[1] + 1])
                                angle_type_lst.append(angle_type)
        m_lst = np.array(molecule_lst)
        molecule_lst = m_lst[m_lst[:, 0].argsort()]

        if len(bond_type_lst) == 0:
            num_bond_types = 1
        else:
            num_bond_types = int(np.max(bond_type_lst))
        if len(angle_type_lst) == 0:
            num_angle_types = 1
        else:
            num_angle_types = int(np.max(angle_type_lst))

        atomtypes = (
            " Start File for LAMMPS \n"
            + "{0:d} atoms".format(len(self._structure))
            + " \n"
            + "{0:d} bonds".format(len(bonds_lst))
            + " \n"
            + "{0:d} angles".format(len(angles_lst))
            + " \n"
            + "{0} atom types".format(len(sorted_species_list))
            + " \n"
            + "{0} bond types".format(num_bond_types)
            + " \n"
            + "{0} angle types".format(num_angle_types)
            + " \n"
        )

        cell_dimensions = self.simulation_cell()

        masses = "Masses" + "\n\n"
        for ic, el_p in enumerate(sorted_species_list):
            mass = self.structure._pse[el_p].AtomicMass
            masses += "{0:3d} {1:f}  # ({2}) \n".format(ic + 1, mass, el_p)

        atoms = "Atoms \n\n"

        # format: atom-ID, molecule-ID, atom_type, q, x, y, z
        format_str = "{0:d} {1:d} {2:d} {3:f} {4:f} {5:f} {6:f}"
        for atom in molecule_lst:
            id_atom, id_mol, id_species = atom
            x, y, z = coords[id_atom]
            ind = np.argwhere(np.array(species_translate_list) == id_species).flatten()[0]
            el_id = self._structure.species[ind].Abbreviation
            atoms += (
                format_str.format(
                    id_atom + 1, id_mol, id_species, q_dict[el_id], x, y, z
                )
                + "\n"
            )

        if len(bonds_lst) > 0:
            bonds_str = "Bonds \n\n"
            for i_bond, id_vec in enumerate(bonds_lst):
                bonds_str += (
                    "{0:d} {1:d} {2:d} {3:d}".format(
                        i_bond + 1, bond_type_lst[i_bond], id_vec[0], id_vec[1]
                    )
                    + "\n"
                )
        else:
            bonds_str = "\n"

        if len(angles_lst) > 0:
            angles_str = "Angles \n\n"
            for i_angle, id_vec in enumerate(angles_lst):
                angles_str += (
                    "{0:d} {1:d} {2:d} {3:d} {4:d}".format(
                        i_angle + 1, angle_type_lst[i_angle], id_vec[0], id_vec[1], id_vec[2]
                    )
                    + "\n"
                )
        else:
            angles_str = "\n"
        return (
            atomtypes
            + "\n"
            + cell_dimensions
            + "\n"
            + masses
            + "\n"
            + atoms
            + "\n"
            + bonds_str
            + "\n"
            + angles_str
            + "\n"
        )

    def structure_charge(self):
        """
        Create atom structure including the atom charges.

        By convention the LAMMPS atom type numbers are chose alphabetically for the chemical species.

        Returns: LAMMPS readable structure.

        """
        atomtypes = (
            "Start File for LAMMPS \n"
            + "{0:d} atoms".format(self._structure.get_number_of_atoms())
            + " \n"
            + "{0} atom types".format(self._structure.get_number_of_species())
            + " \n"
        )

        cell_dimesions = self.simulation_cell()

        masses = "Masses\n\n"

        for ind, obj in enumerate(self._structure.get_species_objects()):
            masses += "{0:3d} {1:f}".format(ind + 1, obj.AtomicMass) + "\n"

        atoms = "Atoms\n\n"

        coords = self.rotate_positions(self._structure)

        el_charge_lst = self._structure.charge
        el_lst = self._structure.get_chemical_symbols()
        el_alphabet_dict = {}
        for ind, el in enumerate(self._structure.get_species_symbols()):
            el_alphabet_dict[el] = ind + 1
        for id_atom, (el, coord) in enumerate(zip(el_lst, coords)):
            id_el = el_alphabet_dict[el]
            dim = self._structure.dimension
            c = np.zeros(3)
            c[:dim] = coord
            atoms += (
                "{0:d} {1:d} {2:f} {3:.15f} {4:.15f} {5:.15f}".format(
                    id_atom + 1, id_el, el_charge_lst[id_atom], c[0], c[1], c[2]
                )
                + "\n"
            )
        return atomtypes + "\n" + cell_dimesions + "\n" + masses + "\n" + atoms + "\n"

    def structure_atomic(self):
        """
        Write routine to create atom structure static file that can be loaded by LAMMPS

        Returns:

        """
        atomtypes = (
            "Start File for LAMMPS \n"
            + "{0:d} atoms".format(len(self._structure))
            + " \n"
            + "{0} atom types".format(len(self._el_eam_lst))
            + " \n"
        )  # '{0} atom types'.format(structure.get_number_of_species()) + ' \n'

        cell_dimesions = self.simulation_cell()

        masses = "Masses\n\n"

        el_struct_lst = self._structure.get_species_symbols()
        el_obj_lst = self._structure.get_species_objects()
        el_dict = {}
        for id_eam, el_eam in enumerate(self._el_eam_lst):
            if el_eam in el_struct_lst:
                id_el = list(el_struct_lst).index(el_eam)
                el = el_obj_lst[id_el]
                el_dict[el] = id_eam + 1
                masses += "{0:3d} {1:f}".format(id_eam + 1, el.AtomicMass) + "\n"
            else:
                # element in EAM file but not used in structure, use dummy for atomic mass
                masses += "{0:3d} {1:f}".format(id_eam + 1, 1.00) + "\n"

        atoms = "Atoms\n\n"

        coords = self.rotate_positions(self._structure)

        el_lst = self._structure.get_chemical_elements()
        for id_atom, (el, coord) in enumerate(zip(el_lst, coords)):
            id_el = el_dict[el]
            dim = self._structure.dimension
            c = np.zeros(3)
            c[:dim] = coord
            atoms += (
                "{0:d} {1:d} {2:.15f} {3:.15f} {4:.15f}".format(
                    id_atom + 1, id_el, c[0], c[1], c[2]
                )
                + "\n"
            )
        return atomtypes + "\n" + cell_dimesions + "\n" + masses + "\n" + atoms + "\n"

    def rotate_positions(self, structure):
        """
        Rotate all atomic positions in given structure according to new Prism cell

        Args:
            structure: Atoms-like object. Should has .positions attribute

        Returns:
            (list): List of rotated coordinates
        """
        prism = UnfoldingPrism(self._structure.cell)
        coords = [prism.pos_to_lammps(position) for position in structure.positions]
        return coords


def write_lammps_datafile(structure, file_name="lammps.data", cwd=None):
    lammps_str = LammpsStructure()
    lammps_str.el_eam_lst = structure.get_species_symbols()
    lammps_str.structure = structure
    lammps_str.write_file(file_name=file_name, cwd=cwd)
