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

pyiron is the integrated development environment for computational materials science. It connects among other things:

* Atomic structure objects – compatible to the `Atomic Simulation Environment (ASE) <https://wiki.fysik.dtu.dk/ase/>`_.

* Atomistic simulation codes – like `LAMMPS <http://lammps.sandia.gov>`_ and `VASP <https://www.vasp.at>`_.

* Feedback Loops – to construct dynamic simulation life cycles.

* Hierarchical data management – interfacing with storage resources like SQL and `HDF5 <https://support.hdfgroup.org/HDF5/>`_.

* Integrated visualization – based on the `NGLview <https://github.com/arose/nglview>`_.

* Interactive simulation protocols - based on `Jupyter notebooks <http://jupyter.org>`_.

* Object oriented job management – for scaling simulation protocols.

pyiron was initially developed in the `Computational Materials Design department <https://www.mpie.de/CM>`_ of `Joerg Neugebauer <https://www.mpie.de/person/43010/2763386>`_ at the `Max Planck Insitut für Eisenforschung (Max Planck Insitute for iron research) <https://www.mpie.de/2281/en>`_ as a framework for ab initio thermo dynamics. In collaboration with the `Interdisciplinary Centre for Advanced Materials Simulation (ICAMS) <http://www.icams.de>`_ the framework was recently extended for high throughput applications resulting in the opensource release of the pyiron.

************
Getting Help
************
Technical issues and bugs should be reported on `Github <https://github.com/pyiron>`_.

***************
Release history
***************

Release 1.0.0 (2018)
====================
* plugin based infrastructure for flexible extension.
* opensource release - licensed under the BSD license.
* installation available on pip and anaconda.
* moved opensource repository to github.

Release 0.9.3 (2017)
====================
* Name changed from PyIron to pyiron
* Fileoperations implemented (move, copy_to and remove).
* NGLview for visualisation.
* Atoms class speedup.
* Serial- and parallelmaster work with the cluster environment.
* S/PHI/nX support added.
* Python 3.6 support added.

Release 0.9.2 (2016)
====================
* Rewirte serial- and parallelmaster.
* Deprecated Qt environment in favor of jupyter.
* Python 3.5 support added.
* Use anaconda as recommended Python environment.
* Switch to Gitlab rather than subversion.

Release 0.9.1 (2015)
====================
* Linux and Mac OS X support added.
* ASE compatible atom and atoms class.

Release 0.0.1 (2011)
====================
* initial version named PyCMW
