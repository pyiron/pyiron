.. _installation_alternatives:

********************************
Alternative Installation Options
********************************
So far we discussed the installation of pyiron on an individual workstation via conda or on a HPC cluster. In the following we focus on developer-specific setups to install pyiron directly from its source. It is recommended to start with a conda installation and then replace only the pyiron version so that conda can still automatically manage all dependencies/environment settings for you. In case this is not possible, e.g. if conda is not allowed on your HPC cluster, then pyiron can be installed directly from the source code.

Install from Source
===================
For development, it is recommended to first create a conda environment containing all of pyiron's dependencies. The dependencies are available in pyiron's `environment.yml <https://github.com/pyiron/pyiron/blob/main/.ci_support/environment.yml>`_ file.

.. code-block:: bash
    git clone https://github.com/pyiron/pyiron.git
    conda env create -f pyiron/environment.yml

If conda is not available on your machine, the next best thing would be to install pyiron and its dependencies via pip.

Using pip
---------
The default installation via pip installs the latest release version of pyiron. So in case your HPC cluster does not support installing pyiron via conda you can install this release version via pip and then continue with the setup of your remote HPC cluster as described above.

.. code-block:: bash

    pip install pyiron

For those who want to test the nightly releases of pyiron which include the latest status of the master branch you can install those via pip as well:

.. code-block:: bash

    pip install --pre pyiron

Using git
---------
To get the latest pyiron version and access changes on development branches pyiron can also be installed via git. For example you can download the pyiron sourcecode to :code:`~/pyiron/software` using:

.. code-block:: bash

    git clone https://github.com/pyiron/pyiron.git ~/pyiron/software

Based on the previous workstation setup your :code:`~/pyiron` directory should contain the following folders:

.. code-block:: bash

   pyiron/
     projects/
     resources/
     software/

To include this version in your :code:`PYTHONPATH` add the following line to your :code:`~/.profile` or :code:`~/.bashrc` configuration:

.. code-block:: bash

    export PYTHONPATH=${HOME}/pyiron/software/:${PYTHONPATH}

When you import pyiron in any python shell or jupyter notebook it should load the version from :code:`~/pyrion/software`. Finally you can switch to other branches using git:

.. code-block:: bash

    git checkout -b main

In this case we switch to the master branch.

Download pyiron Parameter Files
===============================
For source code based installations it is also possible to download the pyiron resources directly from within pyiron. Simply open a python shell and import pyiron:

.. code-block:: python

   > import pyiron
   > pyiron.install()
   >>> It appears that pyiron is not yet configured, do you want to create a default start configuration (recommended: yes). [yes/no]:
   > yes
   > exit()

This command does the following steps in the background:

* Create a :code:`~/.pyiron` config file – with the default settings (for simple installations)

* Create a :code:`~/pyiron/projects` directory – pyiron can only execute calculations within this project directory to prevent any interference with other tools or simulation management solutions.

* Create a :code:`~/pyiron/resources` directory – this directory includes the link to the executables and potentials, sorted by code.

Install Utility
===============
To setup a local lab with pyiron when the internet connection is limited, we provide a classical installer for Windows, macOS X and Linux which is based on the `conda constructor <https://github.com/conda/constructor>`_. If you do not have anaconda installed you can download this installer and get started with just a single `download <https://github.com/pyiron/pyiron-installer/releases>`_.
