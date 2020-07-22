.. _installation:

============
Installation
============

.. note:: **Before you install:** We provide various levels of environments to test pyiron:

          * :ref:`Local Installation (recommeded) <InstallLocal>`: for Windows, Linux or Mac OS X workstation (interface for local VASP executable, support for the latest jupyterlab based GUI)

          * :ref:`Mybinder.org (beta) <InstallBinder>`: test pyiron directly in your browser (no VASP license, only temporary data storage)

          * :ref:`Docker (for demonstration) <InstallDocker>`: requires docker installation (no VASP license, only temporary data storage)


.. _InstallLocal:

*************************************
Workstation Installation (recommeded)
*************************************
Requirements
============
When you start to develop your own simulation protocols we recommend a local installation. Inside the `pyiron anaconda repository <https://anaconda.org/pyiron>`_ we provide precompiled executables for Linux, Mac OS X and Windows with Python 2.7, 3.5, 3.6 and 3.7 and the other packages are available inside the `conda-forge <https://conda-forge.org>`_ channel.

Install pyiron package
======================
As pyiron is written in Python you can install pyiron either via `anaconda <https://www.anaconda.com>`_ (recommended) or via pip.

Install via anaconda (recommended):
------------------------------------------
To install `anaconda <https://www.anaconda.com>`_  you can download the `anaconda distribution <https://www.anaconda.com/download/>`_. Following the installation update to the latest version of conda from `conda-forge <https://conda-forge.org>`_.

.. code-block:: bash

    conda update -c conda-forge conda

After the update of the anaconda environment you can install pyiron using:

.. code-block:: bash

    conda install -c conda-forge pyiron

Install via pip:
----------------
pip is installed on Linux and Mac Os X by default and is included in most Python distributions. To install pyiron via pip type:

.. code-block:: bash

    pip install pyiron

While the anaconda installation already includes the lammps executable, the pip installation requires the user to include a lammps executable named :code:`lmp_serial` for Linux and Mac Os X or :code:`lmp_serial.exe` for windows in their :code:`PATH`.

Visualization
=============
In addition to the pyiron package we recommend installing the `NGLview visualization framework <https://github.com/arose/nglview>`_.

Stable version – for jupyter notebooks (recommended):
-----------------------------------------------------------

.. code-block:: bash

    conda install -c conda-forge nglview
    jupyter nbextension install nglview --py --sys-prefix
    jupyter nbextension enable nglview --py --sys-prefix

Stable version – for jupyter lab
------------------------------------

.. code-block:: bash

    conda install -c conda-forge nodejs nglview
    jupyter labextension install @jupyter-widgets/jupyterlab-manager --no-build
    jupyter labextension install nglview-js-widgets

Simulation code: Lammps
=======================
pyiron supports the simulation codes `VASP <https://www.vasp.at>`_ for DFT calculation and `Lammps <https://lammps.sandia.gov>`_  for molecular dynamics calculation. While VASP requires a separate license and therefore has to be configured by the user, Lammps is available as opensource code and can be installed from anaconda.

For Linux and Mac Os X (for Python 2.7, 3.5, 3.6 and 3.7):
-----------------------------------------------------

.. code-block:: bash

    conda install -c conda-forge lammps

For windows:
------------
.. code-block:: bash

    conda install -c pyiron lammps


Configuration
=============
After the installation of pyiron we need to configure pyiron. The default configuration can be generated automatically. In the terminal, start a new Python session and import pyiron:

.. code-block:: python

   > import pyiron
   > pyiron.install()
   >>> It appears that pyiron is not yet configured, do you want to create a default start configuration (recommended: yes). [yes/no]:
   > yes
   > exit()

The configuration does the following steps in the background:

* Create an :code:`~/.pyiron` config file – with the default settings (for simple installations)

* Create an :code:`~/pyiron/projects` directory – pyiron can only execute calculation within this project directory to prevent any interference, with other tools or simulation management solutions.

* Create an :code:`~/pyiron/resources` directory – this directory includes the link to the executables and potentials, sorted by code. The potentials for lammps are inside :code:`pyiron_lammps` and those for vasp can be placed in :code:`pyiron_vasp`.

First calculation
=================
After the successful configuration you can start your first pyiron calculation. Navigate to the the projects directory and start a jupyter notebook or jupyter lab session correspondingly:

.. code-block:: bash

    cd ~/pyiron/projects
    jupyter notebook

or

.. code-block:: bash

    cd ~/pyiron/projects
    jupyter lab

Open a new jupyter notebook and inside the notebook you can now validate your pyiron calculation by creating a test project, setting up an initial structure of bcc Fe and visualize it using NGLview.

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

Next step
=========
To get a better overview of all the available functionality inside pyiron we recommend the examples provided in the examples section - :doc:`examples`.

.. _IntallCluster:

**********************
Computer Cluster (HPC)
**********************
While the local Installation is designed to scale beyond a single workstation, further multi user extensions are required like:

* `Jupyterhub <https://github.com/jupyterhub/jupyterhub>`_ for managing multiple Jupyter Sessions.

* `PostgreSQL <https://www.postgresql.org>`_ database for scalability.

* Queuing system for job management.

* Access Control lists for sharing files between users.

For further details please open a support request.

.. _InstallBinder:

*******************
Mybinder.org (beta)
*******************

.. warning:: Mybinder.org is currently in beta stage, it should not take longer than a minute to load. We are sorry for the inconvenience.

You can test pyiron on `Mybinder.org (beta) <https://mybinder.org/v2/gh/pyiron/pyiron/master?urlpath=lab>`_, without the need of a local installation. This installation comes with the following limitations:

* No `VASP <https://www.vasp.at>`_ license, DFT calculation can be imported and loaded but the execution is disabled.

* No visualization of atomistic structures using `NGLview <https://github.com/arose/nglview>`_.

* Only temporary data storage, when you leave your session on `Mybinder.org (beta) <https://mybinder.org/v2/gh/pyiron/pyiron/master?urlpath=lab>`_ the environment is reset.

The `Mybinder service <https://mybinder.org>`_ is the most flexible way to test pyiron and get a first impression.
`Start pyiron on MyBinder.org to test your first pyiron examples. <https://mybinder.org/v2/gh/pyiron/pyiron/master?urlpath=lab>`_

.. _InstallDocker:

**************************
Docker (for demonstration)
**************************
Commonly it is easier to install pyiron directly using anaconda following the `Local Installation (Workstation) <InstallLocal>`_ instead of installing Docker. If you already setup Docker on your system, you might still be interested in downloading the pyiron container. While `Mybinder.org (beta) <https://mybinder.org/v2/gh/pyiron/examples.git/master?urlpath=lab>`_ is based on a similar `Docker <https://www.docker.com>`_ image, running the Docker image locally enables more flexibility. In particular the graphical user interface is fully supported in this version. Still the following limitations remain:

* No `VASP <https://www.vasp.at>`_ license, DFT calculation can be imported and loaded but the execution is disabled.

*  Only temporary data storage, when you shutdown your `Docker <https://www.docker.com>`_ instance the environment is reset.

This installation of pyiron is most suitable for presentations. After the local installation of `Docker <https://www.docker.com>`_ there are two versions to choose from stable version based on `jupyter notebooks <http://jupyter.org>`_ and the latest beta version based on `jupyter lab <https://github.com/jupyterlab/jupyterlab>`_. For both versions the first command downloads the image from `Dockerhub <https://hub.docker.com/r/pyiron/pyiron/>`_ and the second command executes it locally.

Docker image with jupyter notebook (stable)
===========================================

.. code-block:: bash

    docker pull pyiron/pyiron:latest

.. code-block:: bash

    docker run -i -t -p 8888:8888 pyiron/pyiron /bin/bash -c "source /srv/conda/envs/notebook/bin/activate; jupyter notebook --notebook-dir=/home/pyiron/ --ip='*' --port=8888"


Docker image with jupyter lab (beta)
====================================

.. code-block:: bash

    docker pull pyiron/pyiron:latest

.. code-block:: bash

    docker run -i -t -p 8888:8888 pyiron/pyiron /bin/bash -c "source /srv/conda/envs/notebook/bin/activate; jupyter lab --notebook-dir=/home/pyiron/ --ip='*' --port=8888"

Connect
=======
After the run command the following line is displayed: Copy/paste this URL into your browser when you connect for the first time, to login with a token:

.. code-block:: bash

    http://localhost:8888/?token=<your_token>

Open the link with your personal jupyter token :code:`<your_token>` in the browser of your choice. Just like the Binder image also the Docker image comes with the examples preinstalled.


.. toctree::
   :maxdepth:2
