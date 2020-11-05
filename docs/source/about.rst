.. _about:

=====
About
=====

************
Introduction
************

.. image:: ../_static/screenshot.png
    :width: 840px
    :height: 525px
    :align: center
    :alt: Screenshot of pyiron running inside jupyterlab.

pyiron is an integrated development environment for implementing, testing, and running simulations in computational materials science. It combines several tools in a common platform:

• Atomic structure objects – compatible to the `Atomic Simulation Environment (ASE) <https://wiki.fysik.dtu.dk/ase/>`_.

• Atomistic simulation codes – like `LAMMPS <http://lammps.sandia.gov>`_ and `VASP <https://www.vasp.at>`_.

• Feedback Loops – to construct dynamic simulation life cycles.

• Hierarchical data management – interfacing with storage resources like SQL and `HDF5 <https://support.hdfgroup.org/HDF5/>`_.

• Integrated visualization – based on `NGLview <https://github.com/arose/nglview>`_.

• Interactive simulation protocols - based on `Jupyter notebooks <http://jupyter.org>`_.

• Object oriented job management – for scaling complex simulation protocols from single jobs to high-throughput simulations.

pyiron (called pyron) is developed in the `Computational Materials Design department <https://www.mpie.de/CM>`_ of `Joerg Neugebauer <https://www.mpie.de/person/43010/2763386>`_ at the `Max Planck Institut für Eisenforschung (Max Planck Institute for iron research) <https://www.mpie.de/2281/en>`_. While its original focus was to provide a framework to develop and run complex simulation protocols as needed for ab initio thermodynamics it quickly evolved into a versatile tool to manage a wide variety of simulation tasks. In 2016 the `Interdisciplinary Centre for Advanced Materials Simulation (ICAMS) <http://www.icams.de>`_ joined the development of the framework with a specific focus on high throughput applications. In 2018 pyiron was released as open-source project.

************
Getting Help
************
Technical issues and bugs should be reported on `Github <https://github.com/pyiron>`_ all other questions can be asked on `stackoverflow using the tag pyiron <https://stackoverflow.com/questions/tagged/pyiron>`_. 

***************
Release history
***************

Release 0.2.0 (2018)
====================
* Implement interactive interface to communicate with codes at runtime.

Release 0.1.0 (2018)
====================
* opensource release - licensed under the BSD license.
* installation available on pip and anaconda.
* moved opensource repository to github.

Release 0.0.9 (2017)
====================
* Name changed from PyIron to pyiron
* Fileoperations implemented (move, copy_to and remove).
* NGLview for visualisation.
* Atoms class speedup.
* Serial- and parallelmaster work with the cluster environment.
* Python 3.6 support added.

Release 0.0.8 (2016)
====================
* Rewirte serial- and parallelmaster.
* Deprecated Qt environment in favor of jupyter.
* Python 3.5 support added.
* Use anaconda as recommended Python environment.
* Switch to Gitlab rather than subversion.

Release 0.0.5 (2015)
====================
* Linux and Mac OS X support added.
* ASE compatible atom and atoms class.

Release 0.0.1 (2011)
====================
* initial version named PyCMW
