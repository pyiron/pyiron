.. _installation_advanced_config:

**********************
Advanced Configuration
**********************
While the conda-based installation is usually sufficient for workstation installations to get started with pyiron, it can be extended to support your own executables, include your own parameter files, support commercial codes like `VASP <https://www.vasp.at>`_ or updating the database performance by switching from `SQLite <https://www.sqlite.org>`_ to `PostgreSQL <https://www.postgresql.org>`_.

Custom Executables and Parameter Files
======================================
pyiron can either be configured using a configuration file named :code:`~/.pyiron` located in the user's home directory or by specifying environment variables. The options are similar either way, so we start with the configuration file. The default configuration file pyiron assumes if it does not find a configuration file is:

.. code-block:: bash

    [DEFAULT]
    PROJECT_CHECK_ENABLED = False
    FILE = ~/pyiron.db
    RESOURCE_PATHS = ${CONDA_PREFIX}/share/pyiron

The first line :code:`[DEFAULT]` defines the current configuration to overwrite the default configuration. The second line :code:`PROJECT_CHECK_ENABLED` disables the project check which enables pyiron to write to the whole file system. The third lines defines the object index to be stored in an SQLite database file :code:`FILE` which is located in the home directory :code:`~/pyiron.db`. It is important to copy the database in case you change the configuration otherwise existing calculation are lost. Finally the :code:`RESOURCE_PATHS` provides the path to the parameter files. Inside pyiron you can check the current configuration using:

.. code-block:: python

    from pyiron_base import Settings
    s = Settings()
    s._configuration

Below, the individual options are explained one by one:

* the :code:`[DEFAULT]` option defines the current :code:`~/.pyiron` configuration to overwrite the default configuration.

* the :code:`RESOURCE_PATHS` option defines the resource path is a list of :code:`;` separated paths where pyiron checks for resource files. A template of such a resource directory is available on `github <https://github.com/pyiron/pyiron-resources>`_ and it can be downloaded as an archive from the `release page <https://github.com/pyiron/pyiron-resources/releases>`_. We recommend to create a folder :code:`~/pyiron/resources` and store the parameter files and links to the executables there. The links are basically shell scripts which can be modified to load modules. By default the conda path is added, therefore there is no need to add it manually.

* the :code:`PROJECT_PATHS` option is similar to the resource path but for storing simulation protocols rather than parameter files. When the :code:`PROJECT_CHECK_ENABLED` option is set to :code:`true` then the read and write access within pyiron is limited to the directories defined in the :code:`PROJECT_PATHS`. Again multiple directories can be separated by :code:`;`. An alternative but outdated name for this option is :code:`TOP_LEVEL_DIRS`.

Besides the general variables in the :code:`~/.pyiron` configuration, the other settings are used to define the database connection. More detailed examples about the configuration can be found below; for now we continue with the configuration of the database. pyiron can use a database to build an index of the HDF5 files on the file system which accelerates job analysis. By default pyiron uses an `SQLite <https://www.sqlite.org>`_ database for this index, but the database can also be disabled or a `PostgreSQL <https://www.postgresql.org>`_ database can be used to improve performance.

* By default the database is defined by the :code:`FILE` option which is equal to the :code:`DATABASE_FILE` option and gives the path to the `SQLite <https://www.sqlite.org>`_ database file. As the `SQLite <https://www.sqlite.org>`_ database is a file-based database, it struggles with parallel access on a shared file system (common for HPC clusters).

* To address this limitation it is possible to disable the database on HPC clusters using the :code:`DISABLE_DATABASE` option by setting it to :code:`true`. This is commonly used when the calculations are only executed on the remote cluster but the analysis is done on a local workstation or a group server which supports an SQL-based database.

* The other database options, namely :code:`TYPE`, :code:`HOST`, :code:`NAME`, :code:`USER`, :code:`PASSWD` and :code:`JOB_TABLE` define the connection details to connect to a PostgreSQL database. Inside pyiron `sqlalchemy <https://www.sqlalchemy.org>`_ is used to support different SQL-based databases, therefore it is also possible to provide the sqlalchemy connection string directly as :code:`CONNECTION`.

* Finally some pyiron installations use a group management component which is currently in development. They might have additional options in their :code:`~/.pyiron` configuration to enable sharing calculations between different users. These options are :code:`VIEWERUSER`, :code:`VIEWERPASSWD` and :code:`VIEWER_TABLE`. As this is a development feature it is not yet fully documented. Basically those are the access details for the global database viewer, which can read the database entries of all users. With this configuration it is possible to load jobs of other users.

In analogy to the :code:`~/.pyiron` configuration file pyiron also supports using environment variables to configure the pyiron installation. The available environment variables are:

* the :code:`PYIRONCONFIG` environment variable defines the location of the :code:`.pyiron` configuration file.

* the :code:`PYIRONRESOURCEPATHS` environment variable defines the :code:`RESOURCE_PATHS` option.

* the :code:`PYIRONPROJECTPATHS` environment variable defines the :code:`PROJECT_PATHS` option.

* the :code:`PYIRONPROJECTCHECKENABLED` environment variable defines the :code:`PROJECT_CHECK_ENABLED` option.

* the :code:`PYIRONDISABLE` environment variable defines the :code:`DISABLE_DATABASE` option.

* the :code:`PYIRONSQLTYPE`, :code:`PYIRONSQLFILE`, :code:`PYIRONSQHOST`, :code:`PYIRONSQLDATABASE`, :code:`PYIRONUSER` and :code:`PYIRONSQLUSERKEY` environment varaibles define the SQL database connection and can also be summarized in the :code:`PYIRONSQLCONNECTIONSTRING` environment variable.

* the :code:`PYIRONSQLVIEWTABLENAME`, :code:`PYIRONSQLVIEWUSER` and :code:`PYIRONSQLVIEWUSERKEY` environment variables define the SQL viewer connection and can also be summarized in the :code:`PYIRONSQLVIEWCONNECTIONSTRING` environment variable.

To further explain the usage of the different parameters, we discuss common use cases in the following:

Use your own Executable for LAMMPS/ S/PHI/nX or GPAW
---------------------------------------------------
To add your own executables or parameter files it is necessary to initialise a user-defined configuration :code:`~/.pyiron`. You can start with a basic configuration like:

.. code-block:: bash

    [DEFAULT]
    FILE = ~/pyiron.db
    PROJECT_PATHS = ~/pyiron/projects
    RESOURCE_PATHS = ~/pyiron/resources

In this case pyiron can only execute calculations in the :code:`~/pyiron/projects` directory. pyiron can’t delete files outside this directory. Next to the projects directory :code:`~/pyiron/projects` we create a resource directory :code:`~/pyiron/resources` to store links to the executables and the corresponding parameter files. Both directories have to be created by the user and in case no :code:`FILE` option is defined pyiron by default creates an `SQLite <https://www.sqlite.org>`_ database in the resource directory. Example resource directories are available on `Github <https://github.com/pyiron/pyiron-resources/tree/main>`_ . Here we just discuss the LAMMPS resource directory as one example.

.. code-block:: bash

    resources/
      lammps/
        bin/
          run_lammps_2020.03.03.sh
          run_lammps_2020.03.03_mpi.sh
        potentials/
          potentials_lammps.csv

The resource directory contains two sub folders :code:`bin` which includes links to the executables and :code:`potentials` which includes links to the interatomic potentials. The links to the executables are shell script which follow the naming convention :code:`run_<code name>_<version>(_<tag>).sh` the :code:`mpi` tag is used to indicate the MPI-enabled executables. If we take a look at the :code:`run_lammps_2020.03.03_mpi.sh` shell script, it contains the following lines:

.. code-block:: bash

    #!/bin/bash
    mpiexec -n $1 --oversubscribe lmp_mpi -in control.inp;

If you are running on a cluster with a module system like `this one <http://modules.sourceforge.net/>`_ and may be a
good idea configure a clean environment that your job can run, e.g.

.. code-block:: bash

   #!/bin/bash
   module purge
   module load lammps/29Oct20
   mpiexec -n $1 --oversubscribe lmp_mpi -in control.inp;

Scripts with the :code:`mpi` tag are called with two parameters the first being the number of cores the second the number of threads, while regular shell scripts do not get any input parameters. By using shell scripts it is easy to link existing executables which might require loading specific modules or setting environment variables. In the same way the parameter files for pyiron are stored in the csv format which makes them human editable. For shared installations we recommend storing the pyiron resources in a shared directory.

Configure VASP
--------------
The `Vienna Ab initio Simulation Package <https://www.vasp.at>`_ is a popular commercial DFT code which is commonly used for large DFT calculations or high-throughput studies. pyiron implements a VASP wrapper but does not provide a VASP license. Therefore users have to compile their own VASP executable and provide their own VASP pseudopotentials (included with the VASP license). An example configuration for VASP in pyiron is available on `Github <https://github.com/pyiron/pyiron-resources/tree/main/vasp>`_:

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

Similar to the LAMMPS resource directory discussed above the VASP resource directory also contains a :code:`bin` diirectory and a :code:`potentials` directory. By adding the :code:`default` tag we can set the default executable, in particular when compiling multiple variants of the same VASP version. Finally the directories :code:`potpaw` and :code:`potpaw_PBE` contain the VASP pseudopotentials, which are included with the VASP license and have to be added by the user.

PostgreSQL Database
===================
To accelerate the pyiron installation it is recommended to use a `PostgreSQL <https://www.postgresql.org>`_ database rather than the default `SQLite <https://www.sqlite.org>`_ database. To configure the database server, the following options can be added to the :code:`~/.pyiron`:

* :code:`TYPE` the typ of the database, while `sqlalchemy <https://www.sqlalchemy.org>`_ supports a wide range of differnet databases `PostgreSQL <https://www.postgresql.org>`_ is recommended and can be selected by setting the type to :code:`Postgres`.

* :code:`HOST` the database host where the database is running.

* :code:`NAME` the name of the database.

* :code:`USER` the database user, in contrast to many other software packages pyiron requires one database user per system user who is using pyiron. The database is only used to store an index of the calculations executed with pyiron, therefore knowledge gained from accessing the database is limited unless the user has also access to the file system.

* :code:`PASSWD` the database user password. While it is a bad practice to store the database password in the configuration file, the database only contains the the job index. Still it is important that the user creates an pyiron specific password and should never store their system user password in the :code:`.pyiron` configuration file.

* :code:`JOB_TABLE` the name of the database table. pyiron is commonly using one table per user.

A typical :code:`.pyiron` configuration with a `PostgreSQL <https://www.postgresql.org>`_ database might look like this:

.. code-block:: bash

    [DEFAULT]
    TYPE = Postgres
    HOST = database.hpc-cluster.university.edu
    NAME = pyiron
    USER = janj
    PASSWD = **********
    JOB_TABLE = jobs_janj
    PROJECT_PATHS = ~/pyiron/projects
    RESOURCE_PATHS = ~/pyiron/resources

Be careful when updating the database configuration as pyiron does not transfer the content of the database automatically.
