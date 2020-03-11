# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
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
        self.__name__ = "VaspMetadyn"
        self.input = MetadynInput()
        self.supported_primitive_constraints = ["bond", "angle", "torsion", "x_pos", "y_pos", "z_pos"]
        self.supported_complex_constraints = ["linear_combination", "norm", "coordination_number"]
        self._constraint_dict = dict()
        self._complex_constraints = dict()
        self.input.incar["LBLUEOUT"] = True

    def set_primitive_constraint(self, name, constraint_type, atom_indices, biased=False, increment=0.0):
        if self.structure is None:
            raise ValueError("The structure has to be set before a dynamic constraint is assigned")
        self._constraint_dict[name] = dict()
        if constraint_type in self.supported_primitive_constraints:
            self._constraint_dict[name]["constraint_type"] = constraint_type
        else:
            raise ValueError("The constraint type '{}' is not supported".format(constraint_type))
        self._constraint_dict[name]["atom_indices"] = atom_indices
        self._constraint_dict[name]["biased"] = biased
        self._constraint_dict[name]["increment"] = increment

    def _set_primitive_constraint_iconst(self, constraint_type, atom_indices, biased=False):
        if self.structure is None:
            raise ValueError("The structure has to be set before a dynamic constraint is assigned")
        constraint_dict = {"bond": "R", "angle": "A", "torsion": "T",
                           "x_pos": "X", "y_pos": "Y", "z_pos": "Z"}
        if constraint_type not in list(constraint_dict.keys()):
            raise ValueError("Use a compatible constraint type")
        line = len(self.input.iconst._dataset["Value"])
        status = 0
        if biased:
            status = 7
        if constraint_type in ["x_pos", "y_pos", "z_pos"]:
            if isinstance(atom_indices, (list, np.ndarray)):
                a_ind = str(self.sorted_indices[atom_indices[0]] + 1)
            else:
                a_ind = str(self.sorted_indices[atom_indices] + 1)
        elif constraint_type in ["bond"]:
            if len(atom_indices) != 2:
                raise ValueError("For this constraint the atom_indices must be a list or numpy array with 2 values")
            a_ind = ' '.join(map(str, self.sorted_indices[atom_indices] + 1))
        elif constraint_type in ["angle"]:
            if len(atom_indices) != 3:
                raise ValueError("For this constraint the atom_indices must be a list or numpy array with 3 values")
            a_ind = ' '.join(map(str, self.sorted_indices[atom_indices] + 1))
        else:
            raise ValueError("The constraint {} is not implemented!".format(constraint_type))
        constraint_string = "{} {} {}".format(constraint_dict[constraint_type], a_ind, status)
        self.input.iconst.set_value(line, constraint_string)

    def set_complex_constraint(self, name, constraint_type, coefficient_dict, biased=False, increment=0.0):
        if self.structure is None:
            raise ValueError("The structure has to be set before a dynamic constraint is assigned")
        if constraint_type not in self.supported_complex_constraints:
            raise ValueError("Use a compatible complex constraint type")
        self._complex_constraints[name] = dict()
        self._complex_constraints[name]["constraint_type"] = constraint_type
        if len(list(coefficient_dict.keys())) != len(list(self._constraint_dict.keys())):
            raise ValueError("The number of coefficients should match the number of primitive contstraints")
        self._complex_constraints[name]["coefficient_dict"] = coefficient_dict
        self._complex_constraints[name]["biased"] = biased
        self._complex_constraints[name]["increment"] = increment

    def _set_complex_constraint(self, constraint_type, coefficients, biased=False):
        constraint_dict = {"linear_combination": "S", "norm": "C", "coordination_number": "D"}
        if constraint_type not in list(constraint_dict.keys()):
            raise ValueError("Use a compatible constraint type")
        status = 0
        if biased:
            status = 5
        coeffs = " ".join(map(str, coefficients))
        constraint_string = "{} {} {}".format(constraint_dict[constraint_type], coeffs, status)
        line = len(self.input.iconst._dataset["Value"])
        self.input.iconst.set_value(line, constraint_string)

    def write_constraints(self):
        increment_list = list()
        linear_constraint_order = list()
        for key, constraint in self._constraint_dict.items():
            linear_constraint_order.append(key)
            self._set_primitive_constraint_iconst(constraint_type=constraint["constraint_type"],
                                                  atom_indices=constraint["atom_indices"], biased=constraint["biased"])
            increment_list.append(constraint["increment"])

        for constraint in self._complex_constraints.values():
            coefficients = [constraint["coefficient_dict"][val] for val in linear_constraint_order]
            self._set_complex_constraint(constraint_type=constraint["constraint_type"],
                                         coefficients=coefficients, biased=constraint["biased"])
            increment_list.append(constraint["increment"])
        self.input.incar["INCREM"] = " ".join(increment_list)

    def write_input(self):
        """
        Call routines that generate the INCAR, POTCAR, KPOINTS and POSCAR input files
        """
        self.write_constraints()
        super(VaspMetadyn, self).write_input()


class MetadynInput(Input):

    def __init__(self):
        super(MetadynInput, self).__init__()
        self.iconst = GenericParameters(input_file_name=None, table_name="iconst", val_only=True, comment_char="!")
        self.penaltypot = GenericParameters(input_file_name=None, table_name="penaltypot", val_only=True,
                                            comment_char="!")

    def write(self, structure, modified_elements, directory=None):
        """
        Writes all the input files to a specified directory

        Args:
            structure (atomistics.structure.atoms.Atoms instance): Structure to be written
            directory (str): The working directory for the VASP run
        """
        # Writing the constraints, increments, and penalty potentials
        super(MetadynInput, self).write(structure, modified_elements, directory)
        self.iconst.write_file(file_name="ICONST", cwd=directory)
        self.penaltypot.write_file(file_name="PENALTYPOT", cwd=directory)

    def to_hdf(self, hdf):
        super(MetadynInput, self).to_hdf(hdf)
        with hdf.open("input") as hdf5_input:
            self.iconst.to_hdf(hdf5_input)
            self.penaltypot.to_hdf(hdf5_input)

    def from_hdf(self, hdf):
        super(MetadynInput, self).from_hdf(hdf)
        with hdf.open("input") as hdf5_input:
            self.iconst.from_hdf(hdf5_input)
            self.penaltypot.from_hdf(hdf5_input)
