# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from collections import OrderedDict
import numpy as np

__author__ = "Sudarsan Surendralal"
__copyright__ = "Copyright 2017, Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Sudarsan Surendralal"
__email__ = "surendralal@mpie.de"
__status__ = "development"
__date__ = "Sep 1, 2017"



class Doscar(object):

    def __init__(self):
        self.doscar_dict = dict()

    def from_file(self, filename="DOSCAR"):

        with open(filename, "r") as f:
            lines = f.readlines()
            num_points = int(lines[5].strip().split()[2])
            print(num_points)
            for line in lines[6: num_points+8]:
                line = line.strip()



    def _parse_total(self):
        pass

