# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
from collections import OrderedDict
import numpy as np
from pyiron.base.generic.parameters import GenericParameters
import decimal as dec

try:
    from ase.calculators.lammps import Prism

except ImportError:
    try:
        from ase.calculators.lammpsrun import Prism
    except ImportError:
        from ase.calculators.lammpsrun import prism as Prism

__author__ = "Joerg Neugebauer, Sudarsan Surendralal, Yury Lysogorskiy, Jan Janssen, Markus Tautschnig"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - " \
                "Computational Materials Design (CM) Department"
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
            super(UnfoldingPrism, self).__init__(cell, pbc=pbc, tolerance=float('1e-{}'.format(digits)))
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
        self.car_prec = dec.Decimal('10.0') ** \
                        int(np.floor(np.log10(max((xhi, yhi, zhi)))) - digits)
        self.dir_prec = dec.Decimal('10.0') ** (-digits)
        self.acc = float(self.car_prec)
        self.eps = np.finfo(xhi).eps

        # For rotating positions from ase to lammps
        apre = np.array(((xhi, 0, 0),
                         (xyp, yhi, 0),
                         (xzp, yzp, zhi)))
        self.R = np.dot(np.linalg.inv(cell), apre)

        # Actual lammps cell may be different from what is used to create R
        eps = 1e-10

        def fold(vec, pvec, i):
            p = pvec[i]
            x = vec[i] + 0.5 * p
            n = (np.mod(x, p) - x) / p
            # print ('prism: ', n, x, p)
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
        
        d_a = apre[0, 0]/2 - apre[1,0]
        if np.abs(d_a) < eps:
            if d_a < 0:
                print ('debug: apply shift')
                apre[1, 0] += 2 * d_a
                apre[2, 0] += 2 * d_a
        self.A = apre
        self.Ainv = np.linalg.inv(self.A)
        self.prism = None

        if self.is_skewed() and \
                (not (pbc[0] and pbc[1] and pbc[2])):
            raise RuntimeError('Skewed lammps cells MUST have '
                               'PBC == True in all directions!')

    def unfold_cell(self, cell):
        """
        Unfold LAMMPS cell to original
        
        Args:
            cell: LAMMPS cell,

        Returns:
            unfolded cell
        """
        a = cell[0]
        bp = cell[1]
        cpp = cell[2]
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
    def __init__(self, input_file_name=None):
        super(LammpsStructure, self).__init__(input_file_name=input_file_name,
                                              table_name="structure_inp",
                                              comment_char="#",
                                              val_only=True)
        self._structure = None
        self._potential = None
        self._el_eam_lst = []
        self.atom_type = None
        self.cutoff_radius = None
        self.digits = 10

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
        if self.atom_type == 'full':
            input_str = self.structure_full()
        elif self.atom_type == 'bond':
            input_str = self.structure_bond()
        elif self.atom_type == 'charge':
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
        input_str = ''
        self.load_string(input_str)

    # def f2s(self, f):
    #     return str(dec.Decimal(repr(f)).quantize(self.car_prec,
    #                                              dec.ROUND_HALF_EVEN))

    def simulation_cell(self):
        """
        
        Returns:

        """
        # dim = self._structure.dimension
        # amat = self._structure.cell

        self.prism = UnfoldingPrism(self._structure.cell, digits=15)
        xhi, yhi, zhi, xy, xz, yz = self.prism.get_lammps_prism_str()
        # Please, be carefull and not round xhi, yhi,..., otherwise you will get too skew cell from LAMMPS.
        # These values are already checked in UnfoldingPrism to fullfill LAMMPS skewness criteria
        simulation_cell = '0. {} xlo xhi\n'.format(xhi) + \
                          '0. {} ylo yhi\n'.format(yhi) + \
                          '0. {} zlo zhi\n'.format(zhi)

        if self.prism.is_skewed():
            simulation_cell += '{0} {1} {2} xy xz yz\n'.format(xy, xz, yz)

        # if s.VERBOSE():
        #     s.warning("triclinic cells not supported for test purposes",
        #               module="lammps.lammps.LammpsStructure.simulation_cell")
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
        # print "bond_type: ", bond_type

        if self.structure.bonds is None:
            if self.cutoff_radius is None:
                bonds_lst = self.structure.get_bonds(max_shells=1)
            else:
                bonds_lst = self.structure.get_bonds(radius=self.cutoff_radius)
            bonds = []
            # id_mol = 0
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

        atomtypes = ' Start File for LAMMPS \n' + \
                    '{0:d} atoms'.format(len(self._structure)) + ' \n' + \
                    '{0:d} bonds'.format(len(bonds)) + ' \n' + \
                    '{0} atom types'.format(self._structure.get_number_of_species()) + ' \n' + \
                    '{0} bond types'.format(np.max(bond_type)) + ' \n'

        cell_dimesions = self.simulation_cell()

        masses = 'Masses \n\n'
        el_obj_list = self._structure.get_species_objects()
        for object_id, el in enumerate(el_obj_list):
            masses += '{0:3d} {1:f}'.format(object_id + 1, el.AtomicMass) + '\n'

        atoms = 'Atoms \n\n'

        # atom_style bond
        # format: atom-ID, molecule-ID, atom_type, x, y, z
        format_str = '{0:d} {1:d} {2:d} {3:f} {4:f} {5:f} '
        if self._structure.dimension == 3:
            for id_atom, (x, y, z) in enumerate(coords):
                id_mol, id_species = 1, el_dict[elements[id_atom]]  # elList[id_atom]
                # print id_atom + 1, id_mol, id_species + 1, x, y, z
                atoms += format_str.format(id_atom + 1, id_mol, id_species + 1, x, y, z) + '\n'
        elif self._structure.dimension == 2:
            for id_atom, (x, y) in enumerate(coords):
                id_mol, id_species = 1, el_dict[elements[id_atom]]  # elList[id_atom]
                # print id_atom + 1, id_mol, id_species + 1, x, y, z
                atoms += format_str.format(id_atom + 1, id_mol, id_species + 1, x, y, 0.) + '\n'
        else:
            raise ValueError("dimension 1 not yet implemented")

        bonds_str = 'Bonds \n\n'
        for i_bond, (i_a, i_b, b_type) in enumerate(bonds):
            bonds_str += '{0:d} {1:d} {2:d} {3:d}'.format(i_bond + 1, b_type, i_a, i_b) + '\n'

        return atomtypes + '\n' + cell_dimesions + '\n' + masses + '\n' + atoms + '\n' + bonds_str + '\n'

    def structure_full(self):
        """
        Write routine to create atom structure static file for atom_type='full' that can be loaded by LAMMPS

        Returns:

        """
        coords = self.rotate_positions(self._structure)

        # extract electric charges from potential file
        q_dict = {}
        for el in self._structure.get_species_symbols():
            q_dict[el] = float(self.potential.get("set group {} charge".format(el)))

        species_translate_list = list()
        sorted_species_list = self._structure.get_species_symbols()
        for el in self._structure.species:
            ind = np.argwhere(sorted_species_list == el.Abbreviation).flatten()[-1]
            species_translate_list.append(ind)

        # analyze structure to get molecule_ids, bonds, angles etc
        molecule_lst, bonds_lst, angles_lst = [], [], []

        # species_lst = structure.get_species_objects()
        # for id_el, el in enumerate(structure.species):
        #     el.id = id
        # el_lst = structure.get_chemical_elements()

        num_atoms_in_molecule = 3
        neighbors = self._structure.get_neighbors(num_neighbors=num_atoms_in_molecule + 2)
        # print "neighbors: ", neighbors.distances
        id_mol = 0
        indices = self._structure.indices
        for id_el, id_species in enumerate(indices):
            el = self._structure.species[id_species]
            # print "id: ", id, el.Abbreviation, neighbors.indices[id][0:2]
            if el.Abbreviation in ["O"]:
                # print "id_mol: ", id_mol
                id_mol += 1
                molecule_lst.append([id_el, id_mol, id_species])
                # Just to ensure that the attached atoms are indeed H atoms
                # id_n1, id_n2 = np.intersect1d(neighbors.indices[id_el], self._structure.select_index("H"))[0:2]
                id_n1, id_n2 = neighbors.indices[id_el][0:2]
                # print "id: ", id, id_n1, len(el_lst), el_lst[1].id
                molecule_lst.append([id_n1, id_mol, species_translate_list[indices[id_n1]]])
                molecule_lst.append([id_n2, id_mol, species_translate_list[indices[id_n2]]])

                bonds_lst.append([id_el + 1, id_n1 + 1])
                bonds_lst.append([id_el + 1, id_n2 + 1])

                angles_lst.append([id_n1 + 1, id_el + 1, id_n2 + 1])
            elif el.Abbreviation not in ["H"]:  # non-bonded ions
                id_mol += 1
                molecule_lst.append([id_el, id_mol, id_species])

        m_lst = np.array(molecule_lst)
        molecule_lst = m_lst[m_lst[:, 0].argsort()]
        # print "m_lst: ", m_lst
        # print "mol: ", molecule_lst

        atomtypes = ' Start File for LAMMPS \n' + \
                    '{0:d} atoms'.format(len(self._structure)) + ' \n' + \
                    '{0:d} bonds'.format(len(bonds_lst)) + ' \n' + \
                    '{0:d} angles'.format(len(angles_lst)) + ' \n' + \
                    '{0} atom types'.format(self._structure.get_number_of_species()) + ' \n' + \
                    '{0} bond types'.format(1) + ' \n' + \
                    '{0} angle types'.format(1) + ' \n'

        cell_dimensions = self.simulation_cell()

        masses = 'Masses' + '\n\n'
        el_obj_list = self._structure.get_species_objects()
        for object_id, el in enumerate(el_obj_list):
            masses += '{0:3d} {1:f}'.format(object_id + 1, el.AtomicMass) + '\n'

        atoms = 'Atoms \n\n'

        # format: atom-ID, molecule-ID, atom_type, q, x, y, z
        format_str = '{0:d} {1:d} {2:d} {3:f} {4:f} {5:f} {6:f}'
        for atom in molecule_lst:
            id_atom, id_mol, id_species = atom
            # print id_atom, id_mol, id_species
            x, y, z = coords[id_atom]
            el_id = self._structure.species[id_species].Abbreviation
            atoms += format_str.format(id_atom + 1, id_mol, id_species + 1, q_dict[el_id], x, y, z) + '\n'

        if len(bonds_lst) > 0:
            bonds_str = 'Bonds \n\n'
            for i_bond, id_vec in enumerate(bonds_lst):
                bonds_str += '{0:d} {1:d} {2:d} {3:d}'.format(i_bond + 1, 1, id_vec[0], id_vec[1]) + '\n'
        else:
            bonds_str = "\n"

        if len(angles_lst) > 0:
            angles_str = 'Angles \n\n'
            for i_angle, id_vec in enumerate(angles_lst):
                # print "id: ", i_angle, id_vec
                angles_str += '{0:d} {1:d} {2:d} {3:d} {4:d}'.format(i_angle + 1, 1, id_vec[0], id_vec[1], id_vec[2]) \
                              + '\n'
        else:
            angles_str = "\n"
        return atomtypes + '\n' + cell_dimensions + '\n' + masses + '\n' + atoms + '\n' \
               + bonds_str + '\n' + angles_str + '\n'

    def structure_charge(self):
        """
        Create atom structure including the atom charges.
        
        By convention the LAMMPS atom type numbers are chose alphabetically for the chemical species.
        
        Returns: LAMMPS readable structure.

        """
        atomtypes = 'Start File for LAMMPS \n' + \
                    '{0:d} atoms'.format(self._structure.get_number_of_atoms()) + ' \n' + \
                    '{0} atom types'.format(self._structure.get_number_of_species()) + ' \n'

        cell_dimesions = self.simulation_cell()

        masses = 'Masses\n\n'
        
        for ind, obj in enumerate(self._structure.get_species_objects()):
            masses += '{0:3d} {1:f}'.format(ind+1, obj.AtomicMass) + '\n'


        atoms = 'Atoms\n\n'

        coords = self.rotate_positions(self._structure)

        el_charge_lst = self._structure.charge
        el_lst = self._structure.get_chemical_symbols()
        el_alphabet_dict = {}
        for ind,el in enumerate(self._structure.get_species_symbols()):
            el_alphabet_dict[el] = ind+1
        for id_atom, (el, coord) in enumerate(zip(el_lst, coords)):
            id_el = el_alphabet_dict[el]
            dim = self._structure.dimension
            c = np.zeros(3)
            c[:dim] = coord
            atoms += '{0:d} {1:d} {2:f} {3:.15f} {4:.15f} {5:.15f}'.format(
                id_atom + 1, id_el, el_charge_lst[id_atom], c[0], c[1], c[2]) + '\n'
        return atomtypes + '\n' + cell_dimesions + '\n' + masses + '\n' + atoms + '\n'

 
    def structure_atomic(self):
        """
        Write routine to create atom structure static file that can be loaded by LAMMPS
        
        Returns:

        """
        atomtypes = 'Start File for LAMMPS \n' + \
                    '{0:d} atoms'.format(len(self._structure)) + ' \n' + \
                    '{0} atom types'.format(len(
                        self._el_eam_lst)) + ' \n'  # '{0} atom types'.format(structure.get_number_of_species()) + ' \n'

        cell_dimesions = self.simulation_cell()

        masses = 'Masses\n\n'

        el_struct_lst = self._structure.get_species_symbols()
        el_obj_lst = self._structure.get_species_objects()
        # el_struct_lst = structure.get_chemical_symbols()
        # el_obj_lst = structure.get_chemical_elements()
        el_dict = {}
        for id_eam, el_eam in enumerate(self._el_eam_lst):
            if el_eam in el_struct_lst:
                id_el = list(el_struct_lst).index(el_eam)
                el = el_obj_lst[id_el]
                el_dict[el] = id_eam + 1
                masses += '{0:3d} {1:f}'.format(id_eam + 1, el.AtomicMass) + '\n'
            else:
                # element in EAM file but not used in structure, use dummy for atomic mass
                masses += '{0:3d} {1:f}'.format(id_eam + 1, 1.00) + '\n'

        atoms = 'Atoms\n\n'

        coords = self.rotate_positions(self._structure)

        el_lst = self._structure.get_chemical_elements()
        for id_atom, (el, coord) in enumerate(zip(el_lst, coords)):
            id_el = el_dict[el]
            dim = self._structure.dimension
            c = np.zeros(3)
            c[:dim] = coord
            atoms += '{0:d} {1:d} {2:.15f} {3:.15f} {4:.15f}'.format(id_atom + 1, id_el, c[0], c[1], c[2]) + '\n'
        return atomtypes + '\n' + cell_dimesions + '\n' + masses + '\n' + atoms + '\n'

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


def write_lammps_datafile(structure, file_name='lammps.data', cwd=None):
    lammps_str = LammpsStructure()
    lammps_str.el_eam_lst = structure.get_species_symbols()
    lammps_str.structure = structure
    lammps_str.write_file(file_name=file_name, cwd=cwd)
