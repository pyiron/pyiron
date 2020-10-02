.. _installation:

============
Installation
============
******************
Conda Installation
******************
The recommended way to install pyiron is via the conda package manager in a Linux environment. So if you are using Windows we recommend using the `Windows subsystem for Linux<https://docs.microsoft.com/en-us/windows/wsl/install-win10>`_ to install pyiron and if you are on Mac Os X we recommend using a `virutal machine/ virtual box<https://www.virtualbox.org>`_. Native installations on both Windows and Mac Os X are possible but are restricted to molecular dynamics calculations with interatomic potentials and do not support density functional theory(DFT) codes. We collaborate with the open-source community at `conda-forge<https://conda-forge.org>`_ to not only provide the pyiron package via their community channel, but also executables for compatible simulation codes like `GPAW<https://wiki.fysik.dtu.dk/gpaw/>`_, `LAMMPS<https://lammps.sandia.gov>`_ and `S/PHI/nX<https://sxrepo.mpie.de>` and their parameter files like pseudo potentials and interatomic potentials. To get started you can install pyiron using: 

.. code-block:: bash

    conda install -c conda-forge pyiron
    

Optional Dependencies 
=====================
All the optional dependencies can also be installed via conda directly to simplify the setup of your simulation environment. 

NGLview for visualising atomistic sturctures
--------------------------------------------
In pyiron we use the `NGLview<http://nglviewer.org/nglview/latest/>`_ package to visualise atomistic structures directly in the jupyter notebook. To enable this feature the NGLview conda package can be installed using: 

.. code-block:: bash

    conda install -c conda-forge nglview

In case you prefer `jupyter lab<https://jupyter.org>`_ over jupyter notebooks, you can also install NGLview for jupyter lab, this requires a few additional dependencies: 

.. code-block:: bash

    conda install -c conda-forge nodejs nglview
    jupyter labextension install @jupyter-widgets/jupyterlab-manager --no-build
    jupyter labextension install nglview-js-widgets

In addition to `NGLview<http://nglviewer.org/nglview/latest/>`_ the first line also installs nodejs which is required to install your own jupyterlab plugins and rebuild jupyter lab. The following two lines install the jupyterlab extensions. Starting with the jupyterlab manager and followed by the NGLview javascript widget. During the installation of `NGLview<http://nglviewer.org/nglview/latest/>`_ it is important to confirm that the NGLview version installed via conda is the same as the version of the NGLview javascript widget: 

.. code-block:: bash

    conda list nglview
    jupyter labextension list

  
LAMMPS for interatomic potentials
---------------------------------
`LAMMPS<https://lammps.sandia.gov>`_ stands for Large-scale Atomic/Molecular Massively Parallel Simulator and it is one of the most popular open-source moelcular dynamics simulation codes for simulating solid-state materials (metals, semiconductors). As part of the pyiron project we maintain the conda package for LAMMPS to simplifiy the installation. Still these executables can be a factor 2-3 slower than executables compiled on the hardware directly. Therefore the conda version is recommended for testing and experimental use your workstation but for high performance computing(HPC) we highly recommend to compile LAMMPS yourself and we discuss how to link your own LAMMPS version below in the addvanced section. 

.. code-block:: bash

    conda install -c conda-forge -c pyiron lammps

On the conda-forge channel we provide the LAMMPS executable for both serial execution and parallel execution using the message passing interface(MPI) as well as python bindings for both serial and parallel execution. The LAMMPS version on the pyiron channel is for native windows installations only and it is limited to serial execution and does not support Python bindings. We therefore highly recommend using the Linux subsystem for Windows rather than the native Windows installation. 

S/PHI/nX for density functional theory
--------------------------------------
The `S/PHI/nX<https://sxrepo.mpie.de>` DFT code is an open-source DFT code developed in close collaboration with the pyiron developers, therefore it is the recommended DFT code to be used with pyiron. The applications of S/PHI/nX range from constrained magnetic calculation to charged defects which makes it suitable for ab initio thermodynamics and beyond. The S/PHI/nX DFT codes is only available on Linux therefore we recommend Windows users to use the Windows subsystem for Linux and Mac Os X users to use a virtual machine. 

.. code-block:: bash

    conda install -c conda-forge sphinxdft

Again the generic executables provided via the conda-forge channel can be a factor of 2-3 slower than compiling S/PHI/nX directly. We therefore recommend to use these executables for testing only. 

GPAW an alternative density functional theory code
--------------------------------------------------
`GPAW<https://wiki.fysik.dtu.dk/gpaw/>`_ is a open-source realspace DFT simulation code implemented in pyiron which is popular because of its Python bindings which allow accessing parameters of the DFT code during the run time. GPAW can be installed on Linux directly via conda:

.. code-block:: bash

    conda install -c conda-forge sphinxdft

Again the generic executables provided via the conda-forge channel can be a factor of 2-3 slower than compiling GPAW directly. We therefore recommend to use these executables for testing only. 

**********************
Advanced configuration
**********************
While the conda based installation is commonly sufficient for workstation installations to get started with pyiron it can be extended to support your own executable, include your own parameter files, support commercial codes like `VASP<https://www.vasp.at>`_ or updating the database performance by switching from `SQLite<https://www.sqlite.org>`_ to `PostgeSQL<https://www.postgresql.org>`_. 

Custom executables and parameter files
======================================
pyiron can either be configured by a configuration file named `~/.pyiron` located in the users home directory or by specifying environment variables. The options are similar so we start with the configuration file. The default configuration file pyiron assums if it does not find a configuration file is:

.. code-block:: bash

    [DEFAULT]  
    PROJECT_CHECK_ENABLED = False
    FILE = ~/pyiron.db
    RESOURCE_PATHS = ${CONDA_PREFIX}/share/pyiron

The first line `[DEFAULT]` defines the current configuration to overwrite the default configuration. The second line `PROJECT_CHECK_ENABLED` disables the project check which enables pyiron to write to the whole file system. The third lines defines the object index to be stored in an SQLite database file `FILE` which is located in the home directory `~/pyiron.db`. It is important to copy the database in case you change the configuration otherwise existing calculation are lost. Finally the `RESOURCE_PATHS` provides the path to the parameter files. Inside pyiron you can check the current configuration using: 

.. code-block:: python

    from pyiron_base import Settings
    s = Settings()
    s._configuration 

In the following the individual options are explained one by one:

* the `[DEFAULT]` option defines the current `~/.pyiron` configuration to overwrite the default configuration.

* the `RESOURCE_PATHS` option defines the resource path is a list of `;` separated paths where pyiron checks for resource files. A template of such a resource directory is available on `github <https://github.com/pyiron/pyiron-resources>`_ and it can be downloaded as an archive from the `release page <https://github.com/pyiron/pyiron-resources/releases>`_. We recommend to create a folder `~/pyiron/resources` and store the parameter files and likes to the executables there. The links are basically shell scripts which can be modified to load modules. By default the conda path is added, therefore there is no need to add it manually. 

* the `PROJECT_PATHS` option is similar to the resource path but for storing simulation protocols rather than parameter files. When the `PROJECT_CHECK_ENABLED` option is set to `true` then the read and write access within pyiron is limited to the directories defined in the `PROJECT_PATHS`. Again multiple directories can be separated by `;`. An alternative but outdated name for this option is `TOP_LEVEL_DIRS`. 

Besides the general variables in the `~/.pyiron` configuration, the other settings are used to define the database connection. More detailed examples about the configuration can be found below, for now we continue with the configuration of the database. pyiron can use a database to build an index of the HDF5 files on the file system which accelerates the analysis. By default pyiron uses an `SQLite<https://www.sqlite.org>`_ database for this index, but the database can also be disabled or a `PostgeSQL<https://www.postgresql.org>`_ database can be used to improve the performance. 

* By default the database is defined by the `FILE` option which is equal to the `DATABASE_FILE` option and gives the path to the `SQLite<https://www.sqlite.org>`_ database file. As the `SQLite<https://www.sqlite.org>`_ database is a file based database it struggles with parallel access on a shared file system like it is typically used in HPC clusters. 

* To address this limitation it is possible to disable the database on HPC clusters using the `DISABLE_DATABASE` option by setting it to `true`. This is commonly used when the calculation are only executed on the remote cluster but the analysis is done on a local workstation or a group server which supports an SQL-based database. 

* The other database options, namely `TYPE`, `HOST`, `NAME`, `USER`, `PASSWD` and `JOB_TABLE` define the connection details to connect to a PostgreSQL database. Inside pyiron `sqlalchemy<https://www.sqlalchemy.org>`_ is used to support different SQL-based databases, therefore it is also possible to provide the sqlalchemy connection string directly as `CONNECTION`. 

* Finally some pyiron installations use a group management component which is currently in development. They might have additional options in their `~/.pyiron` configuration to enable sharing calculations between different users. These options are `VIEWERUSER`, `VIEWERPASSWD` and `VIEWER_TABLE`. As this is a development feature it is not yet fully documented. Basically those are the access details for the global database viewer, which can read the database entries of all users. With this configuration it is possible to load jobs of other users. 

In analogy to the `~/.pyiron` configuration file pyiron also supports using environment variables to configure the pyiron installation. The available environment variables are: 

* the `PYIRONCONFIG` environment variable defines the location of the `.pyiron` configuration file. 

* the `PYIRONRESOURCEPATHS` environment variable defines the `RESOURCE_PATHS` option.

* the `PYIRONPROJECTPATHS` environment variable defines the `PROJECT_PATHS` option.

* the `PYIRONPROJECTCHECKENABLED` environment variable defines the `ROJECT_CHECK_ENABLED` option.

* the `PYIRONDISABLE` environment variable defines the `DISABLE_DATABASE` option.

* the `PYIRONSQLTYPE`, `PYIRONSQLFILE`, `PYIRONSQHOST`, `PYIRONSQLDATABASE`, `PYIRONUSER` and `PYIRONSQLUSERKEY` environment varaibles define the SQL database connection and can also be summarized in the `PYIRONSQLCONNECTIONSTRING` environment variable. 

* the `PYIRONSQLVIEWTABLENAME`, `PYIRONSQLVIEWUSER` and `PYIRONSQLVIEWUSERKEY` environment variables define the SQL viewer connection and can also be summarized in the `PYIRONSQLVIEWCONNECTIONSTRING` environment variable. 

To further explain the usage of the different parameters we discuss common use cases in the following: 

Use you own executable for LAMMPS/ S/PHI/nX or GPAW
---------------------------------------------------
To add your own executables or parameter files it is necessary to initialise a user defined configuration `~/.pyiron`. You can start with a basic configuration like: 

.. code-block:: bash

    [DEFAULT]  
    FILE = ~/pyiron.db
    PROJECT_PATHS = ~/pyiron/projects
    RESOURCE_PATHS = ~/pyiron/resources
    
In this case pyiron can only execute calculations in the `~/pyiron/projects` directory. In particular pyiron users can not delete files outside this directory. Next to the projects directory `~/pyiron/projects` we create a resource directory `~/pyiron/resources` to store links to the executables and the corresponding parameter files. Both directories have to be created by the user and in case no `FILE` option is defined pyiron by default creates an `SQLite<https://www.sqlite.org>`_ database in the resource directory. Example resource directories are available on `Github<https://github.com/pyiron/pyiron-resources/tree/master>`_ . Here we just discuss the LAMMPS resource directory as one example.

.. code-block:: bash

    resources/
      lammps/
        bin/
          run_lammps_2020.03.03.sh
          run_lammps_2020.03.03_mpi.sh
        potentials/
          potentials_lammps.csv

The resource directory contains two sub folders `bin` which includes links to the executables and `potentials` which includes links to the interatomic potentials. The links to the executables are shell script which follow the naming convention `run_<code name>_<version>(_<tag>).sh` the `mpi` tag is used to indicate the MPI-enabled executables. If we take a look at the `run_lammps_2020.03.03_mpi.sh` shell script, it contains the following lines: 

.. code-block:: bash

    #!/bin/bash
    mpiexec -n $1 --oversubscribe lmp_mpi -in control.inp;

Scripts with the `mpi` tag are called with two parameters the first being the number of cores the second the number of threads, while regular shell scripts do not get any input parameters. By using shell scripts it is easy to link exeisting executables which might require loading specific modules or setting environment variables. In the same way the parameter files for pyiron are stored in the csv format which makes them human editable. For shared installations we recommend storing the pyiron resources in a shared directory. 

Configure VASP
--------------
The `Vienna Ab initio Simulation Package<https://www.vasp.at>`_ is a popular commercial DFT code which is commonly used for large DFT calculations or high-throughput studies. pyiron implements a VASP wrapper but does not provide as VASP license. Therefore the users have to compile their own VASP executable and provide their own VASP Pseudopotentials which are included with the VASP license. An example configuration of VASP for pyiron is available on `Github<https://github.com/pyiron/pyiron-resources/tree/master/vasp>`_: 

.. code-block:: bash

    resources/
      vasp/
        bin/
          run_vasp_5.4.4_default.sh
          run_vasp_5.4.4_default_mpi.sh
        potentials/
          potpaw/
          potpaw_PBE/
          potentials_vasp.csv
          potentials_vasp_lda_default.csv
          potentials_vasp_pbe_default.csv

Similar to the LAMMPS resource directory discussed above the VASP resource directory also contains a `bin` diirectory and a `potentials` directory. By adding the `default` tag we can set the default executable, in particular when compiling multiple variants of the same VASP version. Finally the directories `potpaw` and `potpaw_PBE` contain the VASP pseudo potentials, which are included with the VASP license and have to be added by the user. 

PostgreSQL database
===================
To accelerate the pyiron installation it is recommended to use a `PostgeSQL<https://www.postgresql.org>`_ database rather than the default `SQLite<https://www.sqlite.org>`_ database. To configure the database server, the following options can be added to the `~/.pyiron`:

* `TYPE` the typ of the database, while `sqlalchemy<https://www.sqlalchemy.org>`_ supports a wide range of differnet databases `PostgeSQL<https://www.postgresql.org>`_ is recommended and can be selected by setting the type to `Postgres`. 

* `HOST` the database host, where the database is running. 

* `NAME` the name of the database.

* `USER` the database user, in contrast to many other software packages pyiron requires one database user per system user who is using pyiron. The database is only used to store an index of the calculations executed with pyiron, therefore knowledge gained from accessing the database is limited unless the user has also access to the file system. 

* `PASSWD` the database user password. While it is a bad practice to store the database password in the configuration file, the database only contains the the job index. Still it is important that the user creates an pyiron specific password and should never store their system user password in the `.pyiron` configuration file. 

* `JOB_TABLE` the name of the database table. pyiron is commonly using one table per user. 

A typical `.pyiron` configuration with a `PostgeSQL<https://www.postgresql.org>`_ database might look like this: 

.. code-block:: bash

    [DEFAULT]  
    TYPE = Postgres
    HOST = pyiron-db-cmmc.esc.rzg.mpg.de
    NAME = pyiron
    USER = janj
    PASSWD = **********
    JOB_TABLE = jobs_cmmc
    PROJECT_PATHS = ~/pyiron/projects
    RESOURCE_PATHS = ~/pyiron/resources  

Be careful when updating the database configuration as pyiron does not transfer the content of the database automatically. 

Remote HPC cluster
==================
While the previous section discussed the installation of pyiron on a local workstation the following section discusses how to configure a remote HPC cluster to transfer jobs to the HPC cluser for execution and back for analysis. 

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

***************************************
Demonstration and Training environments
***************************************
cloud solutions
===============
You can test pyiron on `Mybinder.org (beta) <https://mybinder.org/v2/gh/pyiron/pyiron/master?urlpath=lab>`_, without the need of a local installation. This installation comes with the following limitations:

* No `VASP <https://www.vasp.at>`_ license, DFT calculation can be imported and loaded but the execution is disabled.

* No visualization of atomistic structures using `NGLview <https://github.com/arose/nglview>`_.

* Only temporary data storage, when you leave your session on `Mybinder.org (beta) <https://mybinder.org/v2/gh/pyiron/pyiron/master?urlpath=lab>`_ the environment is reset.

The `Mybinder service <https://mybinder.org>`_ is the most flexible way to test pyiron and get a first impression.
`Start pyiron on MyBinder.org to test your first pyiron examples. <https://mybinder.org/v2/gh/pyiron/pyiron/master?urlpath=lab>`_

Docker container
================
For demonstration purposes we provide Docker containers on  `Dockerhub <https://hub.docker.com/r/pyiron/pyiron/>`_ these can be downloaded and executed locally once docker is installed. But they have to be used carefully as by default they do not provide any permanent storage, so all information is lost once the docker container is shut down. To download the docker container use: 

.. code-block:: bash

    docker pull pyiron/pyiron:latest

After downloading the docker container you can use it with either jupyter notebook:

.. code-block:: bash

    docker run -i -t -p 8888:8888 pyiron/pyiron /bin/bash -c "source /srv/conda/envs/notebook/bin/activate; jupyter notebook --notebook-dir=/home/pyiron/ --ip='*' --port=8888"

or with jupyter lab:

.. code-block:: bash

    docker run -i -t -p 8888:8888 pyiron/pyiron /bin/bash -c "source /srv/conda/envs/notebook/bin/activate; jupyter lab --notebook-dir=/home/pyiron/ --ip='*' --port=8888"

After the run command the following line is displayed: Copy/paste this URL into your browser when you connect for the first time, to login with a token:

.. code-block:: bash

    http://localhost:8888/?token=<your_token>

Open the link with your personal jupyter token :code:`<your_token>` in the browser of your choice. Just like the Binder image also the Docker image comes with the examples preinstalled.

Install utility
===============
To setup a local lab with pyiron we provide a classical installer for Windows, Mac Os X and Linux which is based on the conda constructor. If you do not have anaconda installed you can download this installer and get started with just a single download. 
https://github.com/pyiron/pyiron-installer/releases

***************
Getting started
***************
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


.. toctree::
   :maxdepth:2
