# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
from phonopy.structure.atoms import PhonopyAtoms
import phonopy.structure.spglib as spg

__author__ = "Jan Janssen"
__copyright__ = "Copyright 2017, Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "development"
__date__ = "Sep 1, 2017"


def analyse_phonopy_equivalent_atoms(atoms):
    """
    Args:
        mode (str): ['total', 'numeric', 'str']
    """

    positions = atoms.positions
    cell = atoms.cell
    types = atoms.get_chemical_symbols()
    types = list(types)
    natom = len(types)
    positions = np.reshape(np.array(positions), (natom, 3))
    cell = np.reshape(np.array(cell), (3, 3))
    for i in range(3):
        xs = np.dot(np.array(positions), np.linalg.inv(np.array(cell)))
        unitcell = PhonopyAtoms(symbols=types, cell=cell, scaled_positions=xs)
        ops = spg.get_symmetry(unitcell)
    return ops['equivalent_atoms']
