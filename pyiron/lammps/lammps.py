# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from pyiron.lammps.interactive import LammpsInteractive

__author__ = "Joerg Neugebauer, Sudarsan Surendralal, Jan Janssen"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - " \
                "- Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Sudarsan Surendralal"
__email__ = "surendralal@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


class Lammps(LammpsInteractive):
    """
    Class to setup and run and analyze LAMMPS simulations which is a derivative of
    atomistics.job.generic.GenericJob. The functions in these modules are written in such the function names and
    attributes are very generic (get_structure(), molecular_dynamics(), version) but the functions are written to handle
    LAMMPS specific input/output.

    Args:
        project (pyiron.project.Project instance):  Specifies the project path among other attributes
        job_name (str): Name of the job

    Attributes:
        input (lammps.Input instance): Instance which handles the input
    """

    def __init__(self, project, job_name):
        super(Lammps, self).__init__(project, job_name)
        self.__name__ = "Lammps"
        self._executable_activate(enforce=True)
