# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from pyiron.base.generic.parameters import GenericParameters
from pyiron.vasp.vasp import Vasp
from pyiron.vasp.base import Input

__author__ = "Sudarsan Surendralal"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Sudarsan Surendralal"
__email__ = "surendralal@mpie.de"
__status__ = "testing"
__date__ = "March 1, 2020"


class VaspMetadyn(Vasp):

    def __init__(self, project, job_name):
        super(VaspMetadyn, self).__init__(project, job_name)
        self.input = MetadynInput()


class MetadynInput(Input):

    def __init__(self):
        super(MetadynInput, self).__init__()
        self.iconst = Iconst()


class Iconst(GenericParameters):
    """
    Class to control the ICONST file of a vasp simulation
    """

    def __init__(self, input_file_name=None, table_name="iconst"):
        super(Iconst, self).__init__(
            input_file_name=input_file_name,
            table_name=table_name,
            val_only=True,
            comment_char="!",
        )

    def load_default(self):
        """
        Loads the default file content
        """
        file_content = """\
ICONST file generated with pyiron
"""
        self.load_string(file_content)


