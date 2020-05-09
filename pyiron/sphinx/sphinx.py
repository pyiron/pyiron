# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import warnings
from pyiron.sphinx.interactive import SphinxInteractive

__author__ = "Osamu Waseda, Jan Janssen"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "development"
__date__ = "Sep 1, 2017"


class Sphinx(SphinxInteractive):
    """
    Class to setup and run Sphinx simulations which is a derivative of pyiron_atomistics.job.generic.GenericJob.
    The functions in these modules are written in such the function names and attributes are very generic
    (get_structure(), molecular_dynamics(), version) but the functions are written to handle Sphinx specific input and
    output.

    Args:
        project: Project object (defines path where job will be created and stored)
        job_name (str): name of the job (must be unique within this project path)
    """

    def __init__(self, project, job_name):
        super(Sphinx, self).__init__(project, job_name)
        self.__name__ = "Sphinx"
        self.__version__ = (
            None
        )  # Reset the version number to the executable is set automatically
        self._executable_activate(enforce=True)
