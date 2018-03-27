.. _installation:

============
Installation
============

The pyiron IDE is designed to scale from single workstations to heterogeneous computing infrastructures. The principles
are the same, still to leverage the scalability of a computer cluster, it is reasonable to setup separate database
servers or login nodes.

---------------------
Windows - workstation
---------------------


install pyiron
==============

from anaconda (recommended)
---------------------------
 ::

    conda install -c pyiron pyiron

from pip
--------
 ::

    pip install pyiron
    
pyiron requires spglib for symmetry analysis of atomistic structures. To compile spglib using pip Microsoft Visual C++
14.0 is required. Please have a look `here <https://wiki.python.org/moin/WindowsCompilers#Which_Microsoft_Visual_C.2B-
.2B-_compiler_to_use_with_a_specific_Python_version_.3F>`_


from source (for developers)
----------------------------

pyiron config (of the .pyiron file)
-----------------------------------

The `.pyiron` file in your home directory must look like this::

    [DEFAULT]
    PROJECT_PATHS = ~/pyiron/projects
    RESOURCE_PATHS = ~/pyiron/resources

-------------------
Linux - workstation
-------------------

install pyiron
==============

from anaconda
-------------
 ::

    conda install -c pyiron pyiron

from pip
--------
 ::

    pip install pyiron
    

from source (for developers)
----------------------------

pyiron config (of the .pyiron file)
-----------------------------------
The .pyiron file in your home directory must look like this::

    [DEFAULT]  
    PROJECT_PATHS = ~/pyiron/projects
    RESOURCE_PATHS = ~/pyiron/resources
    

Simulation codes
----------------

Lammps Library
==============
 ::

   conda install -c pyiron lammps

----------------------
Mac OS X - workstation
----------------------

install pyiron
==============

from anaconda
-------------
 ::

    conda install -c pyiron pyiron

from pip
--------
 ::

   pip install pyiron

from source (for developers)
----------------------------

pyiron config (of the .pyiron file)
-----------------------------------
The .pyiron file in your home directory must look like this::

    [DEFAULT]  
    PROJECT_PATHS = ~/pyiron/projects
    RESOURCE_PATHS = ~/pyiron/resources


---------------
Linux - Cluster
---------------
Based on the Linux workstation installation, we replace the Jupyter notebook server with Jupyterhub and the file-based
SQLite database with Postgres database for scalability. For more details please contact the developers.::

    [DEFAULT]  
    TYPE = Postgres
    HOST = <database host>
    NAME = <database name>
    USER = <username>
    PASSWD = <databasekey>
    JOB_TABLE = jobs_cmmc
    PROJECT_PATHS = /home
    VIEWERUSER = <username read-only user>
    VIEWERPASSWD = <databasekey read-only user>
    VIEWER_TABLE = <database read-only view>
    RESOURCE_PATHS = /pyiron/resources



.. toctree::
   :maxdepth: 1
