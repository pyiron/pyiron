.. _installation:

============
Installation
============


Depending on what you want, there are several options to install and configure pyiron as seen blow. If you first want to quickly see how pyiron feels like, use the binder service on `Mybinder.org (beta) <https://mybinder.org/v2/gh/pyiron/pyiron/main?urlpath=lab>`_ (more info in :ref:`mybinder_explenation`), without the need for a local installation.
For a workstation installation with conda it boils down to :code:`conda install -c conda-forge pyiron`, see :ref:`installation_workstation` for more information and additional packages to be installed.


.. toctree::
   :maxdepth: 2

   installation_workstation
   installation_quickstart
   installation_advanced_config
   installation_hpc_config
   installation_alternatives
   installation_demo_train_envs


Getting Started
***************
Finally once you have installed pyiron you can quickly test your installation with the following minimalistic example. Many more examples are available in the `Github repository <https://github.com/pyiron/pyiron/tree/main/notebooks>`_.

First Calculation
=================
After the successful configuration you can start your first pyiron calculation. Navigate to the the projects directory and start a jupyter notebook or jupyter lab session correspondingly:

.. code-block:: bash

    cd ~/pyiron/projects
    jupyter notebook

or

.. code-block:: bash

    cd ~/pyiron/projects
    jupyter lab

Open a new jupyter notebook and inside the notebook you can now validate your pyiron calculation by creating a test project, setting up an initial structure of bcc Fe, and visualising it using NGLview.

.. code-block:: python

    from pyiron import Project
    pr = Project('test')
    basis = pr.create_structure('Fe', 'bcc', 2.78)
    basis.plot3d()

Finally a first lammps calculation can be executed by:

.. code-block:: python

    ham = pr.create_job(pr.job_type.Lammps, 'lammpstestjob')
    ham.structure = basis
    ham.potential = ham.list_potentials()[0]
    ham.run()

Next Steps
==========
To get a better overview of all the available functionality inside pyiron we recommend the examples provided in the examples section - :doc:`examples`.

.. toctree::
   :maxdepth:2
