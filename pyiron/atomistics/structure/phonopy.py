# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
from phonopy.structure.atoms import PhonopyAtoms
import spglib as spg
from pyiron_base import Settings

__author__ = "Osamu Waseda"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Osamu Waseda"
__email__ = "waseda@mpie.de"
__status__ = "development"
__date__ = "Sep 1, 2018"

s = Settings()


def analyse_phonopy_equivalent_atoms(atoms, symprec=1e-5, angle_tolerance=-1.0):
    """
    Args: (read phonopy.structure.spglib for more details)
        symprec:
            float: Symmetry search tolerance in the unit of length.
        angle_tolerance:
            float: Symmetry search tolerance in the unit of angle deg.
                If the value is negative, an internally optimized routine
                is used to judge symmetry.

    """
    s.publication_add(publication())
    positions = atoms.get_scaled_positions()
    cell = atoms.cell
    types = atoms.get_chemical_symbols()
    types = list(types)
    natom = len(types)
    positions = np.reshape(np.array(positions), (natom, 3))
    cell = np.reshape(np.array(cell), (3, 3))
    unitcell = PhonopyAtoms(symbols=types, cell=cell, scaled_positions=positions)
    ops = spg.get_symmetry(unitcell, symprec=symprec, angle_tolerance=angle_tolerance)
    return ops["equivalent_atoms"]


def publication():
    return {
        "phonopy": {
            "phonopy": {
                "journal": "Scr. Mater.",
                "year": "2015",
                "title": "First principles phonon calculations in materials science",
                "author": ["Togo, A", "Tanaka, I"],
                "pages": "1--5",
                "volume": "108",
                "month": "Nov",
            }
        }
    }
