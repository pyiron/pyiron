.. pyiron documentation master file

Highlights
============
Quick start
.. include:: manual_new/highlights.rst

Introduction
============
pyiron: integrated development environment (IDE) for computational material science

* Concepts
* Simulation lifecycle
* Pyiron objects
.. include:: manual_new/introduction.rst


Installation
============

* Requirements
* Package managers
    * Linux
    * Mac
    * Windows
* Pip
* Git Repository
* .pyiron config
* pyiron_data Folder
* Add executables
* Add cluster configuration
* Test your installation
.. include:: manual_new/installation.rst


Pyiron IDE
============

* Jupyter
* Types
* Loops
* Functions
* Import modules
* numpy
* matplotlib
* pandas
.. include:: manual_new/pyiron_ide.rst


Pyiron objects
============
Project object
* import
* create
* Project path
* list
* open, close, with
* get_item
* remove
* jobtable
* queue table

Structure object
* create structure
* Create surface
* create ase bulk
* Ase structure
* Vacancy (delete)
* Substitutional (indexing)
* Interstitial (add)
* Surface (add vacuum)
* Cutting (numpy)
* Grainboundary (Rotation)
* Neighbors (from bcc to simple cubic)
* Periodic table

Hamilton object
* Lammps
    * set structure
    * select potential
    * modify input (Generic Parameters)
    * minimize
    * molecular dynamics
    * output
    * output files
    * get item, list
    * final structure
    * animation
    * delete a job
    * move a job
    * copy a job
    * non Modal Mode
    * manual mode
    * multi cores
    * submit to the queue
    * check queue status
    * delete job from queue
* VASP
    * import existing VASP calculation
* (Sphinx ) - Part of the initial release?

Metajob object
* Skriptjob
* ListMaster
* SerialMaster -> Convergence test
* ParallelMaster -> Energy Volume curve

.. include:: manual_new/pyiron_objects.rst

Workflows
============

* Convergence test
* Energy Volume curve
* Defect calculation (Vacancy or Interstitial)
* Magnetism
* Monte Carlo
* Band structure
* Pseudo Hydrogen
* Evolutionary Algorithm

.. include:: manual_new/workflows.rst

Publication
============

* How to cite pyiron
* Publication using pyiron

.. include:: manual_new/publication.rst

Multi User environment
============
* Jupyterhub
* view mode
.. include:: manual_new/multiuserinstall.rst

Development
============
How to contribute
* Gitlab repository

Debugging
* Set the debug level
* Report a bug

Develop a code interface
* Hamilton
    * write input
    * collect output
* Meta job
    * run_if_modal
.. include:: manual_new/development.rst

Release notes
============
* Version 1.0.0
.. include:: manual_new/release_notes.rst

Frequently Asked Questions
============
.. include:: manual_new/faq.rst

Code documentation
============
* Module Index

Contact
============
.. include:: manual_new/contact.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
