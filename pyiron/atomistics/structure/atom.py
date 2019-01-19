# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
from pyiron.atomistics.structure.periodic_table import PeriodicTable, ChemicalElement
from pyiron.atomistics.structure.sparse_list import SparseArrayElement
from six import string_types

__author__ = "Joerg Neugebauer, Sudarsan Surendralal"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - " \
                "Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Sudarsan Surendralal"
__email__ = "surendralal@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


def atomproperty(name, doc):
    """Helper function to easily create Atom attribute property."""

    def getter(self):
        return self.get(name)

    def setter(self, value):
        self.set(name, value)

    def deleter(self):
        self.delete(name)

    return property(getter, setter, deleter, doc)


# needed for ASE compatibility
# begin ASE
def abcproperty(index):
    """Helper function to easily create Atom ABC-property."""

    def getter(self):
        spos = self.atoms.get_scaled_positions()
        return spos[self.index][index]

    def setter(self, value):
        spos = self.atoms.get_scaled_positions()
        spos[self.index][index] = value
        self.atoms.set_scaled_positions(spos)

    return property(getter, setter, doc='ABC'[index] + '-coordinate')


def xyzproperty(index):
    """Helper function to easily create Atom XYZ-property."""

    def getter(self):
        return self.position[index]

    def setter(self, value):
        self.position[index] = value

    return property(getter, setter, doc='XYZ'[index] + '-coordinate')
# end ASE


class Atom(SparseArrayElement):
    def __init__(self, symbol='X', position=(0, 0, 0), tag=None, momentum=None, mass=None, magmom=None,
                 charge=None, atoms=None, index=None, pse=None, element=None, **qwargs):
        if element is None and symbol:
            element = symbol
        if tag or momentum or mass or magmom or charge:
            raise ValueError('Not supported parameter used!')
        SparseArrayElement.__init__(self, **qwargs)
        # super(SparseArrayElement, self).__init__(**qwargs)
        # verify that element is given (as string, ChemicalElement object or nucleus number
        if pse is None:
            pse = PeriodicTable()

        if element is None or element=='X':
            if "Z" in qwargs:
                el_symbol = pse.atomic_number_to_abbreviation(qwargs["Z"])
                self._lists['element'] = pse.element(el_symbol)
            else:
                raise ValueError('Need at least element name, Chemical element object or nucleus number')
        else:
            if isinstance(element, string_types):
                el_symbol = element
                self._lists['element'] = pse.element(el_symbol)
            elif isinstance(element, str):
                el_symbol = element
                self._lists['element'] = pse.element(el_symbol)
            elif isinstance(element, ChemicalElement):
                self._lists['element'] = element
            else:
                raise ValueError('Unknown element type')

        self._position = np.array(position)

        # ASE compatibility
        self.index = index
        self._atoms = atoms

    #####  ASE compatibility

    def get(self, name):
        """Get name attribute, return None if not explicitely set."""
        if self._atoms is None:
            return self._position
        return self._atoms.positions[self.index]

    def set(self, name, value):
        """Set name attribute to value."""
        if self._atoms is None:
            self._position = value
        else:
            array = self._atoms.positions
            array[self.index] = value

    def cut_reference_to_atoms(self):
        """Cut reference to atoms object."""
        self._position = self.position
        self.index = None
        self._atoms = None

    @property
    def mass(self):
        return float(self.element.AtomicMass)

    @property
    def symbol(self):
        return self.element.Abbreviation

    @property
    def number(self):
        return self.element.AtomicNumber

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        out_str = "Element: " + self.element.Abbreviation
        for key, val in self._lists.items():
            if key not in ["element"]:
                out_str += ", " + str(key) + ": " + str(val)
        return out_str

    def __eq__(self, other):
        if not (isinstance(other, Atom)):
            return False
        conditions = [np.allclose(self.position, other.position), self.symbol == other.symbol]
        return all(conditions)

    position = atomproperty('position', 'XYZ-coordinates')

    x = xyzproperty(0)
    y = xyzproperty(1)
    z = xyzproperty(2)
