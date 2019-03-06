# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import sys
try:
    from ase.atoms import Atoms as ASEAtoms
except ImportError:
    pass
import pyiron.atomistics.structure.atom as atom
import pyiron.atomistics.structure.atoms as atoms

__author__ = "Joerg Neugebauer"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - " \
                "Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"

sys.modules['ase.atom'] = atom
sys.modules['ase.atoms'] = atoms


def publication():
    return {'ase': """\
@article{ase-paper,
  author={Ask Hjorth Larsen and Jens Jørgen Mortensen and Jakob Blomqvist and Ivano E Castelli and Rune Christensen and Marcin Dułak and Jesper Friis and Michael N Groves and Bjørk Hammer and Cory Hargus and Eric D Hermes and Paul C Jennings and Peter Bjerre Jensen and James Kermode and John R Kitchin and Esben Leonhard Kolsbjerg and Joseph Kubal and Kristen Kaasbjerg and Steen Lysgaard and Jón Bergmann Maronsson and Tristan Maxson and Thomas Olsen and Lars Pastewka and Andrew Peterson and Carsten Rostgaard and Jakob Schiøtz and Ole Schütt and Mikkel Strange and Kristian S Thygesen and Tejs Vegge and Lasse Vilhelmsen and Michael Walter and Zhenhua Zeng and Karsten W Jacobsen},
  title={The atomic simulation environment—a Python library for working with atoms},
  journal={Journal of Physics: Condensed Matter},
  volume={29},
  number={27},
  pages={273002},
  url={http://stacks.iop.org/0953-8984/29/i=27/a=273002},
  year={2017}
}"""}
