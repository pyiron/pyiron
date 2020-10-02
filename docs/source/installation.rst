.. _installation:

============
Installation
============
******************
Conda Installation
******************
The recommended way to install pyiron is via the conda package manager in a Linux environment. So if you are using windows we recommend using the Linux subsystem for Windows to install pyiron and if you are on Mac Os X we recommend using a virutal machine. Native installations on both Windows and Mac Os X are possible but are restricted to molecular dynamics calculations with interatomic potentials and do not support density functional theory(DFT) codes. We collaborate with the open-source community at `conda-forge <https://conda-forge.org>`_ to not only provide the pyiron package via their community channel, but also executables for compatible simulation codes like GPAW, LAMMPS and S/PHI/nX and their parameter files like pseudo potentials and interatomic potentials. To get started you can install pyiron using: 

.. code-block:: bash

    conda install -c conda-forge pyiron
    

Optional Dependencies 
=====================
All the optional dependencies can also be installed via conda directly to simplify the setup of your simulation environment. 

NGLview for visualising atomistic sturctures
--------------------------------------------
In pyiron we use the NGLview package to visualise atomistic structures directly in the jupyter notebook. To enable this feature the NGLview conda package can be installed using: 

.. code-block:: bash

    conda install -c conda-forge nglview

In case you prefer jupyter lab over jupyter notebooks, you can also install NGLview for jupyter lab, this requires a few additional dependencies: 

.. code-block:: bash

    conda install -c conda-forge nodejs nglview
    jupyter labextension install @jupyter-widgets/jupyterlab-manager --no-build
    jupyter labextension install nglview-js-widgets

In addition to NGLview the first line also installs nodejs which is required to install your own jupyterlab plugins and rebuild jupyter lab. The following two lines install the jupyterlab extensions. Starting with the jupyterlab manager and followed by the NGLview javascript widget. During the installation of NGLview it is important to confirm that the NGLview version installed via conda is the same as the version of the NGLview javascript widget: 

.. code-block:: bash

    conda list nglview
    jupyter labextension nglview

  
LAMMPS for interatomic potentials
---------------------------------
LAMMPS stands for Large-scale Atomic/Molecular Massively Parallel Simulator and it is one of the most popular open-source moelcular dynamics simulation codes for simulating solid-state materials (metals, semiconductors). As part of the pyiron project we maintain the conda package for LAMMPS to simplifiy the installation. Still these executables can be a factor 2-3 slower than executables compiled on the hardware directly. Therefore the conda version is recommended for testing and experimental use your workstation but for high performance computing(HPC) we highly recommend to compile LAMMPS yourself and we discuss how to link your own LAMMPS version below in the addvanced section. 

.. code-block:: bash

    conda install -c conda-forge -c pyiron lammps

On the conda-forge channel we provide the LAMMPS executable for both serial execution and parallel execution using the message passing interface(MPI) as well as python bindings for both serial and parallel execution. The LAMMPS version on the pyiron channel is for native windows installations only and it is limited to serial execution and does not support Python bindings. We therefore highly recommend using the Linux subsystem for Windows rather than the native Windows installation. 

S/PHI/nX for density functional theory
--------------------------------------
The S/PHI/nX DFT code is an open-source DFT code developed in close collaboration with the pyiron developers, therefore it is the recommended DFT code to be used with pyiron. The applications of S/PHI/nX range from constrained magnetic calculation to charged defects which makes it suitable for ab initio thermodynamics and beyond. The S/PHI/nX DFT codes is only available on Linux therefore we recommend Windows users to use the Windows subsystem for Linux and Mac Os X users to use a virtual machine. 

.. code-block:: bash

    conda install -c conda-forge sphinxdft

Again the generic executables provided via the conda-forge channel can be a factor of 2-3 slower than compiling S/PHI/nX directly. We therefore recommend to use these executables for testing only. 

GPAW an alternative density functional theory code
--------------------------------------------------
GPAW is a open-source realspace DFT simulation code implemented in pyiron which is popular because of its Python bindings which allow accessing parameters of the DFT code during the run time. GPAW can be installed on Linux directly via conda:

.. code-block:: bash

    conda install -c conda-forge sphinxdft

Again the generic executables provided via the conda-forge channel can be a factor of 2-3 slower than compiling GPAW directly. We therefore recommend to use these executables for testing only. 

**********************
Advanced configuration
**********************
While the conda based installation is commonly sufficient for workstation installations to get started with pyiron it can be extended to support your own executable, include your own parameter files, support commercial codes like VASP or updating the database performance by switching from SQLite to PostgeSQL. 

Custom executables and parameter files
======================================
pyiron can either be configured by a configuration file named `~/.pyiron` located in the users home directory or by specifying environment variables. The configuration file can contain the following options: 

.. code-block:: bash

    [DEFAULT]  
    FILE
    DATABASE_FILE
    TOP_LEVEL_DIRS
    PROJECT_PATHS
    RESOURCE_PATHS
    DISABLE_DATABASE
    PROJECT_CHECK_ENABLED
    TYPE
    HOST
    NAME
    USER
    PASSWD
    JOB_TABLE
    VIEWERUSER
    VIEWERPASSWD
    VIEWER_TABLE
    CONNECTION

The available environment variables are: 

.. code-block:: bash

    PYIRONCONFIG
    PYIRONUSER
    PYIRONRESOURCEPATHS
    PYIRONPROJECTPATHS
    PYIRONSQLCONNECTIONSTRING
    PYIRONSQLTABLENAME
    PYIRONSQLVIEWCONNECTIONSTRING
    PYIRONSQLVIEWTABLENAME
    PYIRONSQLVIEWUSER
    PYIRONSQLVIEWUSERKEY
    PYIRONSQLFILE
    PYIRONSQHOST
    PYIRONSQLTYPE
    PYIRONSQLUSERKEY
    PYIRONSQLDATABASE
    PYIRONPROJECTCHECKENABLED
    PYIRONDISABLE

PostgreSQL database
===================
As pyiron provides direct access to a python shell and the python shell provides direct access to the file system we rely on the file system permissions to secure the user directories. 

Remote HPC cluster
==================

.. code-block:: bash

    cluster_primary: cmti001
    cluster:
      cmti001: cmti001.yaml
      cmti002: cmti002.yaml

.. code-block:: bash
    queue_type: REMOTE
    queue_primary: cm
    ssh_host: cmti001.bc.rzg.mpg.de
    ssh_username: janj
    known_hosts: ~/.ssh/known_hosts
    ssh_key: ~/.ssh/id_rsa
    ssh_remote_config_dir: /u/system/SLES12/soft/pyiron/dev/pyiron-resources-cmmc/queues/
    ssh_remote_path: /cmmc/u/janj/remote/
    ssh_local_path: /Users/jan/pyiron/projects/
    ssh_continous_connection: True
    queues:
      cm: {cores_max: 1200, cores_min: 1, run_time_max: 5760}
      cmfe: {cores_max: 1200, cores_min: 1, run_time_max: 5760}
      cmti: {cores_max: 1200, cores_min: 1, run_time_max: 5760}
      s_cmfe: {cores_max: 40, cores_min: 1, run_time_max: 5760}
      p_cmfe: {cores_max: 1200, cores_min: 40, run_time_max: 5760}
      
.. code-block:: bash

    #!/bin/bash
    {%- if cores < 40 %}
    #SBATCH --partition=s.cmfe
    {%- else %}
    #SBATCH --partition=p.cmfe
    {%- endif %}
    #SBATCH --ntasks={{cores}}
    #SBATCH --constraint='[swe1|swe1|swe2|swe3|swe4|swe5|swe6|swe7]'
    {%- if run_time_max %}
    #SBATCH --time={{run_time_max // 60}}
    {%- endif %}
    {%- if memory_max %}
    #SBATCH --mem={{memory_max}}
    {%- else %}
    {%- if cores < 40 %}
    #SBATCH --mem-per-cpu=3GB
    {%- endif %}
    {%- endif %}
    #SBATCH --output=time.out
    #SBATCH --error=error.out
    #SBATCH --job-name={{job_name}}
    #SBATCH --chdir={{working_directory}}
    #SBATCH --get-user-env=L
    {{command}}
   
********************************   
Alternative installation options
********************************
install from source 
===================
using pip
---------
.. code-block:: bash

    pip install pyiron

.. code-block:: bash

    pip install --pre pyiron

using git
---------
.. code-block:: bash

    git clone https://github.com/pyiron/pyiron.git

.. code-block:: bash

    git checkout -b master

setup pyiron configuration
==========================
Again create your `~/.pyiron` configuration file 

download pyiron parameter files
===============================
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

Demonstration and Training environments
=======================================
cloud solutions
---------------
You can test pyiron on `Mybinder.org (beta) <https://mybinder.org/v2/gh/pyiron/pyiron/master?urlpath=lab>`_, without the need of a local installation. This installation comes with the following limitations:

* No `VASP <https://www.vasp.at>`_ license, DFT calculation can be imported and loaded but the execution is disabled.

* No visualization of atomistic structures using `NGLview <https://github.com/arose/nglview>`_.

* Only temporary data storage, when you leave your session on `Mybinder.org (beta) <https://mybinder.org/v2/gh/pyiron/pyiron/master?urlpath=lab>`_ the environment is reset.

The `Mybinder service <https://mybinder.org>`_ is the most flexible way to test pyiron and get a first impression.
`Start pyiron on MyBinder.org to test your first pyiron examples. <https://mybinder.org/v2/gh/pyiron/pyiron/master?urlpath=lab>`_

Docker container
----------------
For demonstration purposes we provide Docker containers on  `Dockerhub <https://hub.docker.com/r/pyiron/pyiron/>`_ these can be downloaded and executed locally once docker is installed. But they have to be used carefully as by default they do not provide any permanent storage, so all information is lost once the docker container is shut down. To download the docker container use: 

.. code-block:: bash

    docker pull pyiron/pyiron:latest

After downloading the docker container you can use it with either jupyter notebook:

.. code-block:: bash

    docker run -i -t -p 8888:8888 pyiron/pyiron /bin/bash -c "source /srv/conda/envs/notebook/bin/activate; jupyter notebook --notebook-dir=/home/pyiron/ --ip='*' --port=8888"

or with jupyter lab:

.. code-block:: bash

    docker run -i -t -p 8888:8888 pyiron/pyiron /bin/bash -c "source /srv/conda/envs/notebook/bin/activate; jupyter lab --notebook-dir=/home/pyiron/ --ip='*' --port=8888"

Install utility
---------------
To setup a local lab with pyiron we provide a classical installer for Windows, Mac Os X and Linux which is based on the conda constructor. If you do not have anaconda installed you can download this installer and get started with just a single download. 
https://github.com/pyiron/pyiron-installer/releases

*****************
OLD DOCUMENTATION
*****************

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
