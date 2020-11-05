# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import sys

try:
    from ase.atoms import Atoms as ASEAtoms
except ImportError:
    pass

__author__ = "Joerg Neugebauer"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


def publication():
    return {
        "ase": {
            "ase-paper": {
                "author": [
                    "Ask Hjorth Larsen",
                    "Jens Jørgen Mortensen",
                    "Jakob Blomqvist",
                    "Ivano E Castelli",
                    "Rune Christensen",
                    "Marcin Dułak",
                    "Jesper Friis",
                    "Michael N Groves",
                    "Bjørk Hammer",
                    "Cory Hargus",
                    "Eric D Hermes",
                    "Paul C Jennings",
                    "Peter Bjerre Jensen",
                    "James Kermode",
                    "John R Kitchin",
                    "Esben Leonhard Kolsbjerg",
                    "Joseph Kubal",
                    "Kristen Kaasbjerg",
                    "Steen Lysgaard",
                    "Jón Bergmann Maronsson",
                    "Tristan Maxson",
                    "Thomas Olsen",
                    "Lars Pastewka",
                    "Andrew Peterson",
                    "Carsten Rostgaard",
                    "Jakob Schiøtz",
                    "Ole Schütt",
                    "Mikkel Strange",
                    "Kristian S Thygesen",
                    "Tejs Vegge",
                    "Lasse Vilhelmsen",
                    "Michael Walter",
                    "Zhenhua Zeng",
                    "Karsten W Jacobsen",
                ],
                "title": "The atomic simulation environment—a Python library for working with atoms",
                "journal": "Journal of Physics: Condensed Matter",
                "volume": "29",
                "number": "27",
                "pages": "273002",
                "url": "http://stacks.iop.org/0953-8984/29/i=27/a=273002",
                "year": "2017",
            }
        }
    }
