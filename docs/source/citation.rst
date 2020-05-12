.. _citation:

======
Citing
======
The pyiron integrated development environment (IDE) for computational materials science - pyiron IDE - is based on a flexible plugin infrastructure. So depending on which modules are used please cite the corresponding papers.

*****************************
pyiron paper (accepted)
*****************************

.. code-block:: bibtex

  @article{pyiron-paper,
    title = {pyiron: An integrated development environment for computational materials science},
    journal = {Computational Materials Science},
    volume = {163},
    pages = {24 - 36},
    year = {2019},
    issn = {0927-0256},
    doi = {https://doi.org/10.1016/j.commatsci.2018.07.043},
    url = {http://www.sciencedirect.com/science/article/pii/S0927025618304786},
    author = {Jan Janssen and Sudarsan Surendralal and Yury Lysogorskiy and Mira Todorova and Tilmann Hickel and Ralf Drautz and Jörg Neugebauer},
    keywords = {Modelling workflow, Integrated development environment, Complex simulation protocols},
  }

For all the other modules/ plugins in particular those hosted at https://gitlab.mpcdf.mpg.de/pyiron (MPIE internal) please ask the developers for the corrsponding references. We try to publish those under the open source license when the initial papers are published. Afterwards they are going to be added to the official `Github repository <https://github.com/pyiron>`_.

**************
external paper
**************
Some of the features in pyiron rely on external codes which should be cited separatly. In alphabetical order:

ASE
===
pyiron is compatible with the `Atomic Simulation Environment (ASE) <https://wiki.fysik.dtu.dk/ase/index.html>`_ structure classes, allowing the user to generate structures using the `ASE framework <https://wiki.fysik.dtu.dk/ase/index.html>`_ and run the simulation within pyiron.

.. code-block:: bibtex

  @article{ase-paper,
    author={Ask Hjorth Larsen and Jens Jørgen Mortensen and Jakob Blomqvist and Ivano E Castelli and Rune Christensen and Marcin Dułak and Jesper Friis and Michael N Groves and Bjørk Hammer and Cory Hargus and Eric D Hermes and Paul C Jennings and Peter Bjerre Jensen and James Kermode and John R Kitchin and Esben Leonhard Kolsbjerg and Joseph Kubal and Kristen Kaasbjerg and Steen Lysgaard and Jón Bergmann Maronsson and Tristan Maxson and Thomas Olsen and Lars Pastewka and Andrew Peterson and Carsten Rostgaard and Jakob Schiøtz and Ole Schütt and Mikkel Strange and Kristian S Thygesen and Tejs Vegge and Lasse Vilhelmsen and Michael Walter and Zhenhua Zeng and Karsten W Jacobsen},
    title={The atomic simulation environment—a Python library for working with atoms},
    journal={Journal of Physics: Condensed Matter},
    volume={29},
    number={27},
    pages={273002},
    url={http://stacks.iop.org/0953-8984/29/i=27/a=273002},
    year={2017}
  }

LAMMPS
======
The `LAMMPS molecular dynamics simulator <http://lammps.sandia.gov>`_ is the default molecular dynamics code used by pyiron.

.. code-block:: bibtex

  @article{lammps,
    title = {Fast Parallel Algorithms for Short-Range Molecular Dynamics},
    journal = {Journal of Computational Physics},
    volume = {117},
    number = {1},
    pages = {1-19},
    year = {1995},
    issn = {0021-9991},
    doi = {https://doi.org/10.1006/jcph.1995.1039},
    url = {http://www.sciencedirect.com/science/article/pii/S002199918571039X},
    author = {Steve Plimpton}
  }

VASP
====
The `Vienna Ab initio Simulation Package <https://www.vasp.at>`_ is the default ab initio used by pyiron.

.. code-block:: bibtex

  @article{Kresse1993,
    title = {Ab initio molecular dynamics for liquid metals},
    author = {Kresse, G. and Hafner, J.},
    journal = {Phys. Rev. B},
    volume = {47},
    issue = {1},
    pages = {558--561},
    numpages = {0},
    month = {Jan},
    publisher = {American Physical Society},
    doi = {10.1103/PhysRevB.47.558},
    url = {https://link.aps.org/doi/10.1103/PhysRevB.47.558}
  }

.. code-block:: bibtex

  @article{Kresse1996a,
    title = {Efficiency of ab-initio total energy calculations for metals and semiconductors using a plane-wave basis set},
    journal = {Computational Materials Science},
    volume = {6},
    number = {1},
    pages = {15-50},
    year = {1996},
    issn = {0927-0256},
    doi = {https://doi.org/10.1016/0927-0256(96)00008-0},
    url = {http://www.sciencedirect.com/science/article/pii/0927025696000080},
    author = {Kresse, G. and Furthm\"uller, J.}
  }

.. code-block:: bibtex

  @article{Kresse1996b,
    title = {Efficient iterative schemes for ab initio total-energy calculations using a plane-wave basis set},
    author = {Kresse, G. and Furthm\"uller, J.},
    journal = {Phys. Rev. B},
    volume = {54},
    issue = {16},
    pages = {11169--11186},
    numpages = {0},
    year = {1996},
    month = {Oct},
    publisher = {American Physical Society},
    doi = {10.1103/PhysRevB.54.11169},
    url = {https://link.aps.org/doi/10.1103/PhysRevB.54.11169}
  }

.. toctree::
   :maxdepth:2
