.. _installation:

============
Installation
============
pyiron is build and tested for Python 2.7, 3.5 and 3.6. It can be installed as a basic python package but we highly recommend installing jupyter notebooks and nglview to have a more modern interface with integrated visualisation capabilities. 

----------------------
Docker (demonstration)
----------------------

The easiest way to test pyiron is to use the docker image. Docker is a framework for operating-system-level virtualization, allowing users to test new software in separated virtual environments. To read more about Docker, visit https://docs.docker.com. After you installed docker on your system, you can download and start the pyiron image using the following commands:  
 ::

    docker pull pyiron/pyiron
    docker run -i -t -p 8888:8888 pyiron/pyiron /bin/bash -c "/opt/conda/bin/jupyter lab --notebook-dir=/opt/notebooks --ip='*' --port=8888 --no-browser --allow-root"

With the first line docker downloads the pyiron image from the pyiron repository https://hub.docker.com/r/pyiron/pyiron/ with the second line docker starts an instance of jupyterlab running pyiron. Ending with the following statement: 
 ::

   Copy/paste this URL into your browser when you connect for the first time,
       to login with a token:
           http://localhost:8888/?token=eb1394c21574a59249b6d36eab4484f59b7f13516f23f152

Copy the line http://localhost:8888/?token=eb1394c21574a59249b6d36eab4484f59b7f13516f23f152 from your terminal to the web browser to access jupyterlab. 



------------------
Local installation
------------------

The local installation of pyiron is designed for a single user running pyiron on a workstation. The installation is operation system independent. We recommend the Anaconda Python distribution https://www.anaconda.com but pyiron can also be installed using pip. For anaconda users, install pyiron: 
 ::

    conda install -c pyiron -c conda-forge nglview=0.6.2.3 jupyter_contrib_nbextensions=0.3.3 pyiron

Alternatively install pyiron from pip (for spglib which is a pyiron dependency a local compiler is required). 
 ::

    pip install nglview==0.6.2.3 jupyter_contrib_nbextensions==0.3.3 pyiron

To validate pyiron is successfully installed, open a Python shell and execute: 
 ::

    import pyiron

The import creates a '~/.pyiron' configuration file and the folders '~/pyiron/projects' and '~/pyiron/resources'. All pyrion projects should be started in the '~/pyiron/projects' folder. pyiron is tracking the pyiron objects within this folder. The '~/pyiron/resources' includes the resources for the individual pyiron plugins. A basic template for the resource directory is avialable at https://github.com/pyiron/pyiron-resources which you can download as a zip file https://github.com/pyiron/pyiron-resources/archive/master.zip . Copy the folders 'pyrion_atomistics', 'pyiron_lammps' and 'pyiron_vasp' to '~/pyiron/resources'.

Afterwards we can test the pyiron visualisation by opening a terminal, navigating to the pyiron projects folder and start a jupyter notebook session: 
 ::
    
    cd ~/pyiron/projects
    jupyter notebook 
    
The jupyter navigator should be started automatically inside the browser, so you can create a new jupyter notebook and enter the following lines: 
 ::
     
    from pyiron import Project
    pr = Project('test')
    basis = pr.create_structure('Fe', 'bcc', 2.78)
    basis.plot3d() 
    
The code above creates a two atom iron bcc structure with a lattice constant of 2.78 and visualizes the structure. To execute a first pyiron calculation we need to add an interface to the simulation code. For lammps this can be done by editing the file '~/pyiron/resources/pyiron_lammps/bin/run_lammps_<version number>.bat' for windows or '~/pyiron/resources/pyiron_lammps/bin/run_lammps_<version number>.sh' for linux / MacOs. The version number is used as an identifier to support multiple versions of the same executable. Sample scripts are provided in https://github.com/pyiron/pyiron-resources . 

In addition to the executables additional resources like emperical potentials '~/pyiron/resources/pyiron_lammps/potentials/' can be stored for each individual code in their resource directory. 

After the executable is configured the first calculation can be executed, using the atomistic structure from above we run: 
 ::
     
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
