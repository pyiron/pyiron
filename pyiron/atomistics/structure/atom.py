# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
from pyiron.atomistics.structure.periodic_table import PeriodicTable, ChemicalElement
from pyiron.atomistics.structure.sparse_list import SparseArrayElement
from six import string_types
from ase.atom import Atom as ASEAtom

__author__ = "Sudarsan Surendralal"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Sudarsan Surendralal"
__email__ = "surendralal@mpie.de"
__status__ = "production"
__date__ = "Aug 1, 2020"


class Atom(ASEAtom, SparseArrayElement):
    """
    Class for representing a single atom derived from the `ASE atom class`_.

    Args:
        symbol (str/pyiron.atomistics.structure.periodic_table.ChemcicalElement): Symbol or elecment object
        position (list/numpy.ndarray): Position of atom in cartesian coordinates
        tag (str): Tag assigned to structure
        momentum (float): Momentum
        mass (float): Atomic mass in a.u.
        magmom (float): Magnetic moment in Bohn Magneton
        charge (float): Charge in e
        atoms (ase.atoms.Atoms): Assigned atoms
        index (int): Assigned index

    .. _ASE atom class: https://wiki.fysik.dtu.dk/ase/ase/atom.html

    """
    def __init__(
        self,
        symbol="X",
        position=(0, 0, 0),
        tag=None,
        momentum=None,
        mass=None,
        magmom=None,
        charge=None,
        atoms=None,
        index=None,
        pse=None,
        element=None,
        **qwargs
    ):
        if element is None:
            element = symbol

        SparseArrayElement.__init__(self, **qwargs)
        # super(SparseArrayElement, self).__init__(**qwargs)
        # verify that element is given (as string, ChemicalElement object or nucleus number
        if pse is None:
            pse = PeriodicTable()

        if element is None or element == "X":
            if "Z" in qwargs:
                el_symbol = pse.atomic_number_to_abbreviation(qwargs["Z"])
                self._lists["element"] = pse.element(el_symbol)
            else:
                raise ValueError(
                    "Need at least element name, Chemical element object or nucleus number"
                )
        else:
            if isinstance(element, string_types):
                el_symbol = element
                self._lists["element"] = pse.element(el_symbol)
            elif isinstance(element, str):
                el_symbol = element
                self._lists["element"] = pse.element(el_symbol)
            elif isinstance(element, ChemicalElement):
                self._lists["element"] = element
            else:
                raise ValueError("Unknown element type")

        # KeyError handling required for user defined elements
        try:
            ASEAtom.__init__(
                self,
                symbol=symbol,
                position=position,
                tag=tag,
                momentum=momentum,
                mass=mass,
                magmom=magmom,
                charge=charge,
                atoms=atoms,
                index=index)
        except KeyError:
            symbol = pse.Parent[symbol]
            ASEAtom.__init__(
                self,
                symbol=symbol,
                position=position,
                tag=tag,
                momentum=momentum,
                mass=mass,
                magmom=magmom,
                charge=charge,
                atoms=atoms,
                index=index)

        # ASE compatibility for tags
        for key, val in qwargs.items():
            self.data[key] = val

    @property
    def mass(self):
        """
        Gives the atomic mass of the atom

        Returns:
            float: The atomic mass in a.u.

        """
        return float(self.element.AtomicMass)

    @property
    def symbol(self):
        """
        The chemical symbol of the atom

        Returns:
            str: The chemical symbol of the atom

        """
        return self.element.Abbreviation

    @property
    def number(self):
        """
        The atomic number of the atom

        Returns:
            int: The atomic number according to the periodic table

        """
        return self.element.AtomicNumber

    def __eq__(self, other):
        if not (isinstance(other, Atom)):
            return False
        conditions = [
            np.allclose(self.position, other.position),
            self.symbol == other.symbol,
        ]
        return all(conditions)


def ase_to_pyiron(ase_obj):
    """
    Convert an ase.atom.Atom object to its equivalent pyiron structure

    Args:
        ase_obj(ase.atom.Atom): The ase atoms instance to convert

    Returns:
        pyiron.atomistics.structure.atom.Atom: The equivalent pyiron Atom

    """
    return Atom(symbol=ase_obj.symbol,
                position=ase_obj.position,
                tag=ase_obj.tag,
                momentum=ase_obj.momentum,
                mass=ase_obj.mass,
                magmom=ase_obj.magmom,
                charge=ase_obj.charge,
                atoms=ase_obj.atoms,
                index=ase_obj.index)
