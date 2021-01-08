# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from collections import OrderedDict
import numpy as np
import scipy.constants
from pyiron_atomistic.atomistics.structure.atoms import Atoms
from pyiron_atomistic.atomistics.structure.periodic_table import PeriodicTable

__author__ = "Sudarsan Surendralal, Jan Janssen"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Sudarsan Surendralal"
__email__ = "surendralal@mpie.de"
__status__ = "production"
__date__ = "Feb 4, 2018"

BOHR_TO_ANGSTROM = (
    scipy.constants.physical_constants["Bohr radius"][0] / scipy.constants.angstrom
)


def read_atoms(filename="structure.sx"):
    """
    Args:
        filename (str): Filename of the sphinx structure file

    Returns:
        pyiron_atomistic.objects.structure.atoms.Atoms instance

    """
    file_string = []
    with open(filename) as f:
        for line in f:
            line = line.strip()
            file_string.append(line)
    cell_trigger = "cell"
    cell_string = list()
    species_list = list()
    species_trigger = "element"
    positions_dict = OrderedDict()
    positions = list()
    pse = PeriodicTable()
    for i, line in enumerate(file_string):
        if cell_trigger in line:
            for j in range(len(file_string)):
                line_str = file_string[i + j]
                cell_string.append(line_str)
                if ";" in line_str:
                    break
        if species_trigger in line:
            species = (
                line.strip().split("=")[-1].replace(";", "").replace('"', "").strip()
            )
            species_list.append(pse.element(species))
            positions_dict[species] = 0
            for j in range(len(file_string) - i):
                line_str = file_string[i + j]
                k = 0
                if "atom" in line_str:
                    break_loop = False
                    while not break_loop:
                        position_string = " ".join(
                            file_string[i + j + k].split("=")[-1]
                        )
                        replace_list = ["[", "]", ";", "}",
                            "movable", "X", "Y", "Z"]
                        for rep in replace_list:
                            position_string = (
                                "".join(position_string).replace(rep, " ").split()
                            )
                        positions.append(
                            np.array(position_string[0].split(","), dtype=float)
                        )
                        positions_dict[species] += 1
                        k += 1
                        if (i + j + k) <= len(file_string) - 1:
                            if (
                                "element" in file_string[i + j + k]
                                or "atom" not in file_string[i + j + k]
                            ):
                                break_loop = True
                    break
    indices = list()
    for i, val in enumerate(positions_dict.values()):
        indices.append(np.ones(val, dtype=int) * i)
    indices = np.hstack(indices)
    replace_list = ["cell", "=", "[", "]", ",", ";"]
    for rep in replace_list:
        cell_string = " ".join(cell_string).replace(rep, " ").split()
    cell = np.array(cell_string, dtype=float).reshape((3, 3)) * BOHR_TO_ANGSTROM
    atoms = Atoms(
        species=species_list,
        indices=indices,
        cell=cell,
        positions=np.array(positions) * BOHR_TO_ANGSTROM,
    )
    return atoms
