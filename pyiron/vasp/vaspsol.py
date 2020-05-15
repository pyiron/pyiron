# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from pyiron.vasp.vasp import Vasp

__author__ = "Jan Janssen"
__copyright__ = "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - " \
                "Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "development"
__date__ = "Sep 1, 2018"


class VaspSol(Vasp):
    def __init__(self, project, job_name):
        super(VaspSol, self).__init__(project, job_name)
        self.__name__ = "VaspSol"
        self.__version__ = None  # Reset the version number to the executable is set automatically
        self._executable = None
        self._executable_activate()
