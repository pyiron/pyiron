# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
from scipy.integrate import cumtrapz

__author__ = "Sudarsan Surendralal"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Sudarsan Surendralal"
__email__ = "surendralal@mpie.de"
__status__ = "production"
__date__ = "Mar 1, 2020"


class Report(object):
    """
    This module is used to parse VASP REPORT files
    """

    def __init__(self):
        self.parse_dict = dict()

    def from_file(self, filename="REPORT"):
        """
        Reads values from files and stores it in the `parse_dict` attribute

        Args:
            filename (str): Path to the file that needs to be parsed
        """
        with open(filename, "r") as f:
            lines = f.readlines()
        rel_lines = [lines[i + 2] for i, line in enumerate(lines) if "Blue_moon" in line]
        if len(rel_lines) > 0:
            [lam, _, _, _] = [val for val in np.genfromtxt(rel_lines, usecols=[1, 2, 3, 4]).T]
            rel_lines = [lines[i] for i, line in enumerate(lines) if "cc>" in line]
            cv = np.genfromtxt(rel_lines, usecols=[2])
            fe = cumtrapz(lam, cv)
            self.parse_dict["cv_full"] = cv
            self.parse_dict["derivative"] = lam
            self.parse_dict["cv"] = cv[:-1]
            self.parse_dict["free_energy"] = fe
