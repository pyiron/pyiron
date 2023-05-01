.. _installation:

============
Install pyiron so you can submit to remote HPCs from a local machine (laptop/workstation, Linux/WSL OS)
============

1. If you have already installed pyiron on your cluster, and it works, we can proceed. If not, click [here](https://pyiron.readthedocs.io/en/latest/source/installation_quickstart.html) and finish that first.

2. To install pyiron on your local machine , first install :code-block: `mamba` via:

.. code-block:: bash 

    wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh

3. Execute the script to install mambaforge:

.. code-block:: bash

    bash Mambaforge-Linux-x86_64.sh

4. Now install pyiron via:

.. code-block:: bash

    bash Mambaforge-Linux-x86_64.sh

5. Now you have the option of installing pyiron in an environment with:

.. code-block:: bash

    mamba create -n YOURENVNAME

Change :code:`YOURENVNAME` to your liking.

6. Then activate your environment with:

.. code-block:: bash

    mamba activate YOURENVNAME

7. Call this to install pyiron:

.. code-block:: bash

    mamba install -c conda-forge pyiron

