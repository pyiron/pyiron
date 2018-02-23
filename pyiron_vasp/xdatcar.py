# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np

"""
This module is used to parse VASP XDATCAR files.
"""

__author__ = "Sudarsan Surendralal"
__copyright__ = "Copyright 2017, Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Sudarsan Surendralal"
__email__ = "Sudarsan Surendralal"
__status__ = "development"
__date__ = "Sep 1, 2017"


class Xdatcar(object):

    """
    This module is used to parse VASP XDATCAR files. This gives the unwrapped coordinates unlike vasprun.xml which wraps
    the coordinates based on the supercell dimensions. This is important in analyzing MD trajectories where
    """

    def __init__(self):
        self.atoms = None