.. _installation:

============
Installation
============

*******************
Before you install:
*******************
We provide various levels of environments to test pyiron: 

* Mybinder.org (beta): test pyiron directly in your browser (no VASP license, no visualization, only temporary data storage) 

* Docker Image: for demonstration (no VASP license, only temporary data storage)

* Anaconda/ pip installation: for Windows, Linux or Mac OS X workstation (interface for local VASP executable, support for the latest jupyterlab based GUI) 

************************
Mybinder.org (first try)
************************
You can test pyiron on Mybinder.org (beta), without the need of a local installation. This installation comes with the following limitations: 

* No VASP license, DFT calculation can be imported and loaded but the execution is disabled.

* No visualization of atomistic structures using NGLview. 

* Only temporary data storage, when you leave your session on Mybinder (beta) the environment is reset. 

The Mybinder service is the most flexible way to test pyiron and get a first impression. 
Start pyiron on MyBinder.org to test your first pyiron examples <button> .

**************************
Docker (for demonstration)
**************************
While Mybinder.org is based on a similar Docker image, running the docker image locally enables more flexibility. In particular the graphical user interface is fully supported in this version. Still the following limitations remain: 

* No VASP license, DFT calculation can be imported and loaded but the execution is disabled. 

*  Only temporary data storage, when you shutdown your Docker instance the environment is reset.

This installation of pyiron is most suitable for presentations. After the local installation of Docker there are tow versions to choose from stable version based on jupyter notebooks and the latest beta version based on jupyter lab. For both versions the first command downloads the image from Dockerhub and the second command executes it locally. 

Docker image with jupyter notebook (stable)
===========================================

.. code-block:: bash

    docker pull pyiron/pyiron:notebook
    docker run -i -t -p 8888:8888 pyiron/pyiron /bin/bash -c "/opt/conda/bin/jupyter notebook --notebook-dir=/opt/notebooks --ip='*' --port=8888 --no-browser --allow-root"

Docker image with jupyter lab (beta)
====================================

.. code-block:: bash

    docker pull pyiron/pyiron:latest
    docker run -i -t -p 8888:8888 pyiron/pyiron /bin/bash -c "/opt/conda/bin/jupyter lab --notebook-dir=/opt/notebooks --ip='*' --port=8888 --no-browser --allow-root"

Connect
=======
After the run command the following line is displayed: 
Copy/paste this URL into your browser when you connect for the first time, to login with a token:

.. code-block:: bash 

    http://localhost:8888/?token=<your_token>

Open the link with your personal jupyter token `<your_token>` in the browser of your choice. Just like the Binder image also the Docker image comes with the examples preinstalled. 

********************************
Local Installation (Workstation)
********************************
Requirements
============
When you start to develop your own simulation protocols we recommend switching to a local installation. For the local installation pyiron requires: 

* Python 2.7+

* Pytables 

* h5py 

* spglib 

* SQLite 

* ASE

Inside the pyiron anaconda repository we provide precompiled executables for Linux, Mac OS X and Windows with Python 2.7, 3.5 and 3.6 and the other packages are available inside the conda-forge channel. 

Install pyiron package 
======================
When your system fulfills these requirements you can install pyiron either via anaconda (recommended) or via pip. 

Install pyiron via anaconda (recommended): 
------------------------------------------

.. code-block:: bash

    conda install -c pyiron -c conda-forge pyiron lammps

Install via pip: 
----------------

.. code-block:: bash 

    pip install pyiron

While the anaconda installation already includes the lammps executable, the pip installation requires the user to include a lammps executable named lmp_serial or lmp_serial.exe for windows in their PATH. 

Visualization 
=============
In addition to the pyiron package we recommend installing the NGLview visualization framework. 

Stable version 0.6.3 – for jupyter notebooks (recommended):
-----------------------------------------------------------

.. code-block:: bash 

    conda install -c conda-forge nglview=0.6.2.3 jupyter_contrib_nbextensions=0.3.3

Test version 1.1.2 – for jupyter lab
------------------------------------

.. code-block:: bash 
    
    conda install ipywidgets=7.1.2 nodejs -c conda-forge
    pip install nglview==1.1.2
    nglview enable
    conda install jupyterlab=0.31.12 -y -c conda-forge
    jupyter-labextension install @jupyter-widgets/jupyterlab-manager@0.33.2
    jupyter-labextension install nglview-js-widgets@1.1.2

Configuration
=============
After the installation of pyiron we need to configure pyiron. The default configuration can be generated automatically. In the terminal, start a new Python session and import pyiron: 
python

.. code-block:: python

   import pyiron
   >>> No pyiron installation found, should pyiron be installed [yes/no]:
   yes
   exit()
   
The configuration does the following steps in the background: 

* Create an ~/.pyiron config file – with the default settings (for simple installations)

* Create an ~/pyiron/projects directory – pyiron can only execute calculation within this project directory to prevent any interference, with other tools or simulation management solutions. 

* Create an ~/pyiron/resources directory – this directory includes the link to the executables and potentials, sorted by code. The potentials for lammps are inside pyiron_lammps and those for vasp can be placed in pyiron_vasp. 

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

.. code-block:: bash

    ham = pr.create_job(pr.job_type.Lammps, 'lammpstestjob')
    ham.structure = basis
    ham.potential = ham.list_potentials()[0]
    ham.run()

Next step
=========
To get a better overview of all the available functionality inside pyiron we recommend the examples provided in the examples section. 

----------------------
Computer Cluster (HPC)
----------------------
While the local Installation is designed to scale beyond a single workstation, further multi user extensions are required like: 

* Jupyterhub for managing multiple Jupyter Sessions.

* PostgreSQL database for scalability.

* Queuing system for job management.

* Access Control lists for sharing files between users.

For further details please open a support request.  


============
Installation
============

.. |pyironbinder| image:: https://mybinder.org/badge.svg
    :target: https://mybinder.org/v2/gh/pyiron/pyiron-docker.git/binder

.. seealso:: Before you install pyiron locally test it live on MyBinder.org |pyironbinder|  

pyiron is build and tested for Python 2.7, 3.5 and 3.6. It can be installed as a basic python package but we highly recommend installing jupyter notebooks and nglview to have a more modern interface with integrated visualisation capabilities. 


----------------------
Docker (demonstration)
----------------------

The easiest way to test pyiron is to use the docker image. Docker is a framework for operating-system-level virtualization, allowing users to test new software in separated virtual environments. To read more about Docker, visit https://docs.docker.com. After you installed docker on your system, you can download and start the pyiron image using the following commands:  

.. code-block:: bash
   
    docker pull pyiron/pyiron
    docker run -i -t -p 8888:8888 pyiron/pyiron /bin/bash -c "/opt/conda/bin/jupyter lab --notebook-dir=/opt/notebooks --ip='*' --port=8888 --no-browser --allow-root"

With the first line docker downloads the pyiron image from the pyiron repository https://hub.docker.com/r/pyiron/pyiron/ with the second line docker starts an instance of jupyterlab running pyiron. Ending with the following statement: 

.. code-block:: bash
   
    Copy/paste this URL into your browser when you connect for the first time, to login with a token:
    http://localhost:8888/?token=<your_token>

Copy the line http://localhost:8888/?token=<your_token> from your terminal to the web browser to access jupyterlab.



------------------
Local installation
------------------

The local installation of pyiron is designed for a single user running pyiron on a workstation. The installation is operation system independent. We recommend the Anaconda Python distribution https://www.anaconda.com but pyiron can also be installed using pip. For anaconda users, install pyiron: 

.. code-block:: bash
   
    conda install -c pyiron -c conda-forge nglview=0.6.2.3 jupyter_contrib_nbextensions=0.3.3 pyiron

Alternatively install pyiron from pip (for spglib which is a pyiron dependency a local compiler is required). 

.. code-block:: bash
   
    pip install nglview==0.6.2.3 jupyter_contrib_nbextensions==0.3.3 pyiron

To validate pyiron is successfully installed, open a Python shell and execute: 

.. code-block:: python
   
    import pyiron
    >>> No pyiron installation found, should pyiron be installed [yes/no]:

If you answer with 'yes' pyiron creates a :code:`~/.pyiron` configuration file and the folders :code:`~/pyiron/projects` and :code:`~/pyiron/resources`. All pyrion projects should be started in the :code:`~/pyiron/projects` folder. pyiron is tracking the pyiron objects within this folder. The :code:`~/pyiron/resources` includes the resources for the individual pyiron plugins. A basic template for the resource directory is avialable at https://github.com/pyiron/pyiron-resources which is downloaded automatically and the folders :code:`pyrion_atomistics`, :code:`pyiron_lammps` and :code:`pyiron_vasp` are copied to :code:`~/pyiron/resources`.

Afterwards we can test the pyiron visualisation by opening a terminal, navigating to the pyiron projects folder and start a jupyter notebook session: 

.. code-block:: bash
   
    cd ~/pyiron/projects
    jupyter notebook 
    
The jupyter navigator should be started automatically inside the browser, so you can create a new jupyter notebook and enter the following lines: 

.. code-block:: python
   
    from pyiron import Project
    pr = Project('test')
    basis = pr.create_structure('Fe', 'bcc', 2.78)
    basis.plot3d() 

If NGLview is not loaded correctly or the atomistic structure is not displayed, which some times happens in particular with the miniconda distribution, it might be reasonable to try the latest development version: 

.. code-block:: bash

    conda install ipywidgets=7.1.2 nodejs -c conda-forge
    pip install nglview==1.1.2
    nglview enable 
    conda install jupyterlab=0.31.12 -y -c conda-forge
    jupyter-labextension install @jupyter-widgets/jupyterlab-manager@0.33.2
    jupyter-labextension install nglview-js-widgets@1.1.2
 
The code above creates a two atom iron bcc structure with a lattice constant of 2.78 and visualizes the structure. To execute a first pyiron calculation we need to add an interface to the simulation code. For lammps this can be done by editing the file :code:`~/pyiron/resources/pyiron_lammps/bin/run_lammps_<version number>.bat` for windows or :code:`~/pyiron/resources/pyiron_lammps/bin/run_lammps_<version number>.sh` for linux / MacOs. The version number is used as an identifier to support multiple versions of the same executable. Sample scripts are provided in https://github.com/pyiron/pyiron-resources . 

In addition to the executables additional resources like emperical potentials :code:`~/pyiron/resources/pyiron_lammps/potentials/` can be stored for each individual code in their resource directory. 

To install lammps http://lammps.sandia.gov on Linux or MacOs we provide a serial executable via anaconda: 

.. code-block:: bash

    conda install lammps -c pyiron
    
The windows executable is available at http://rpm.lammps.org/windows.html . After the installation it has to be linked in :code:`~/pyiron/resources/pyiron_lammps/bin/run_lammps_<version number>.bat`.  

After the executable is configured the first calculation can be executed, using the atomistic structure from above we run: 

.. code-block:: python
   
    ham = pr.create_job(pr.job_type.Lammps, 'lammpstestjob')
    ham.structure = basis 
    ham.potential = ham.list_potentials()[0]
    ham.run()

Congratulation, you executed the first pyiron calculation on your system. 

--------------------
Cluster installation
--------------------

In contrast to the single user installation the cluster installation is using a PostgreSQL database instead of a SQLite database and jupyterhub instead of just jupyter. More details comming soon. 



.. toctree::
   :maxdepth:2
