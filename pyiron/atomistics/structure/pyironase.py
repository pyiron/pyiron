# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import sys
try:
    from ase.atoms import Atoms as ASEAtoms
except ImportError:
    pass
import pyiron.atomistics.structure.atom as atom
import pyiron.atomistics.structure.atoms as atoms

__author__ = "Joerg Neugebauer"
__copyright__ = "Copyright 2017, Max-Planck-Institut für Eisenforschung GmbH - " \
                "Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"

sys.modules['ase.atom'] = atom
sys.modules['ase.atoms'] = atoms
