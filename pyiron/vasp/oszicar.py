# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.


__author__ = "Sudarsan Surendralal"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Sudarsan Surendralal"
__email__ = "surendralal@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2020"


class Oszicar(object):
    """
    This module is used to parse VASP OSZICAR files.

    Attributes:

        parse_dict (dict): A dictionary with all the useful quantities parsed from an OSZICAR file after from_file() is
                           executed

    """

    def __init__(self):
        self.parse_dict = dict()

    def from_file(self, filename="OSZICAR"):
        with open(filename, "r") as f:
            lines = f.readlines()

    def get_energy_pot(self):
        pass
