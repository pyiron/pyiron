# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from collections import OrderedDict

import numpy as np

from pyiron.dft.waves.electronic import ElectronicStructure

__author__ = "Sudarsan Surendralal"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - " \
                "Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Sudarsan Surendralal"
__email__ = "surendralal@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


class Procar(object):
    """
    This module contains routines to parse VASP PROCAR files.
    """

    def __init__(self):
        self._is_spin_polarized = False
        self.dos_dict = OrderedDict()

    def from_file(self, filename):
        with open(filename, "r") as f:
            es_obj = ElectronicStructure()
            lines = f.readlines()
            details_trigger = "# of k-points:"
            details_ready = False
            kpoint_trigger = "k-point"
            band_trigger = "band"
            num_atoms = 0
            for i, line in enumerate(lines):
                line = line.strip()
                if details_trigger in line:
                    num_kpts, num_bands, num_atoms = self._get_details(line)
                    details_ready = True
                if details_ready:
                    if kpoint_trigger in line.split():
                        kpt, weight = self._get_kpoint_details(line)
                        es_obj.add_kpoint(value=kpt, weight=weight)

                    if band_trigger in line.split():
                        eigenvalue, occupancy = self._get_band_details(line)
                        es_obj.kpoints[-1].add_band(eigenvalue=eigenvalue, occupancy=occupancy)
                        band_obj = es_obj.kpoints[-1].bands[-1]
                        band_obj.resolved_dos_matrix, band_obj.orbital_resolved_dos, band_obj.atom_resolved_dos = \
                            self._get_dos_matrix(lines[i + 2: i + num_atoms + 4])
        return es_obj

    @staticmethod
    def _check_if_spin_polarized(line):
        pass

    @staticmethod
    def _get_details(line):
        lst = line.split()
        num_kpts = int(lst[3])
        num_bands = int(lst[7])
        num_atoms = int(lst[11])
        return num_kpts, num_bands, num_atoms

    @staticmethod
    def _get_kpoint_details(line):
        line = line.replace("-", " -")
        lst = line.split()
        kpt = [float(lst[i]) for i in range(4, 7)]
        weight = float(lst[9])
        return kpt, weight

    @staticmethod
    def _get_band_details(line):
        lst = line.split()
        eigval = float(lst[4])
        occ = float(lst[7])
        return eigval, occ

    @staticmethod
    def _get_dos_matrix(lines):
        num_orbitals = len((lines[0].strip()).split()) - 2
        num_atoms = len(lines) - 2
        dos_matrix = np.zeros((num_atoms, num_orbitals))
        orbital_resolved_dos = list()
        atom_resolved_dos = list()
        count = 0
        for i, line in enumerate(lines):
            line = line.strip()
            lst = line.split()

            if i not in [0, len(lines) - 1]:
                dos_matrix[count, :] = np.array([float(val) for val in lst[1: len(lst) - 1]])
                count += 1
                atom_resolved_dos.append(float(lst[-1]))
            if i == len(lines) - 1:
                orbital_resolved_dos = [float(val) for val in lst[1: len(lst) - 1]]

        atom_resolved_dos = np.array(atom_resolved_dos)
        orbital_resolved_dos = np.array(orbital_resolved_dos)
        return dos_matrix, orbital_resolved_dos, atom_resolved_dos
