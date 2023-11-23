.. _installation_workstation:

=============================
Installation on a workstation
=============================

The recommended way to install pyiron is via the conda package manager in a Linux environment. So if you are using Windows we recommend installing the `Windows subsystem for Linux <https://docs.microsoft.com/en-us/windows/wsl/install-win10>`_ before you install pyiron and if you are on macOS X we recommend using a `virtual machine/ virtual box <https://www.virtualbox.org>`_. Native installations on both Windows and macOS X are possible but are restricted to molecular dynamics calculations with interatomic potentials and do not support density functional theory(DFT) codes. We collaborate with the open-source community at `conda-forge <https://conda-forge.org>`_ to not only provide the pyiron package via their community channel, but also executables for compatible simulation codes like `GPAW <https://wiki.fysik.dtu.dk/gpaw/>`_, `LAMMPS <https://lammps.sandia.gov>`_ and `S/PHI/nX <https://sxrepo.mpie.de>`_ and their parameter files like pseudopotentials and interatomic potentials. To get started you can install pyiron using:

.. code-block:: bash

    conda install -c conda-forge pyiron


Optional Dependencies
=====================
All the optional dependencies can also be installed via conda directly to simplify the setup of your simulation environment.

NGLview (Atomistic Structure Visualisation)
--------------------------------------------
In pyiron we use the `NGLview <http://nglviewer.org/nglview/latest/>`_ package to visualise atomistic structures directly in the jupyter notebook. To enable this feature, install NGLview:

.. code-block:: bash

    conda install -c conda-forge nglview

In case you prefer `jupyter lab <https://jupyter.org>`_ over jupyter notebooks, you can also install NGLview for jupyter lab. In older versions this requires a few additional dependencies:

.. code-block:: bash

    conda install -c conda-forge nodejs nglview
    jupyter labextension install @jupyter-widgets/jupyterlab-manager --no-build
    jupyter labextension install nglview-js-widgets

In addition to `NGLview <http://nglviewer.org/nglview/latest/>`_ the first line also installs nodejs which is required to install your own jupyterlab plugins and rebuild jupyter lab. The following two lines install the jupyterlab extensions. Starting with the jupyterlab manager and followed by the NGLview javascript widget. During the installation of `NGLview <http://nglviewer.org/nglview/latest/>`_ it is important to confirm that the NGLview version installed via conda is the same as the version of the NGLview javascript widget:

.. code-block:: bash

    conda list nglview
    jupyter labextension list

Supported simulation packages (quantum engines)
-----------------------------------------------

The following packages are supported to work out-of-the-box with pyiron,
but must be installed independently either using conda or manual compilation. Manually compiled executables can be as much as 2-3x faster than conda-installed executables, and are therefore *strongly* recommended for high performance computing (HPC) usage. We discuss how to link any "homemade" executables to your pyiron installation in the advanced section.

LAMMPS (Molecular Dynamics with Interatomic Potentials)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
`LAMMPS <https://lammps.sandia.gov>`_ stands for Large-scale Atomic/Molecular Massively Parallel Simulator and it is one of the most popular open-source molecular dynamics simulation codes for simulating solid-state materials (metals, semiconductors). As part of the pyiron project we maintain the conda package for LAMMPS to simplifiy its installation.

.. code-block:: bash

    # serial + parallel, for linux and mac systems
    conda install -c conda-forge lammps

    # only serial (no python bindings), for native windows
    conda install -c conda-forge -c pyiron lammps

On the conda-forge channel we provide LAMMPS executables for both serial and parallel (MPI) execution as well as their respective python bindings. The LAMMPS version on the pyiron channel is for native windows installations only and it is limited to serial execution with no Python bindings. We therefore highly recommend using the Linux subsystem for Windows rather than the native Windows installation.

S/PHI/nX (Density Functional Theory)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The `S/PHI/nX <https://sxrepo.mpie.de>`_ DFT code is an open-source DFT code developed in close collaboration with the pyiron developers, therefore it is the recommended DFT code to be used with pyiron. The applications of S/PHI/nX range from constrained magnetic calculations to charged defects which makes it suitable for ab initio thermodynamics and beyond. The S/PHI/nX DFT code is only officially supported for Linux, so we recommend the use of a Linux subsystem (on Windows) or a virtual machine (on mac).

.. code-block:: bash

    conda install -c conda-forge sphinxdft

GPAW (Density Functional Theory)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
pyiron also supports `GPAW <https://wiki.fysik.dtu.dk/gpaw/>`_, an open-source realspace DFT simulation code which is popular because of its Python bindings which allow accessing parameters of the DFT code during the run time. GPAW can be installed on Linux directly via conda:

.. code-block:: bash

    conda install -c conda-forge gpaw


Additional simulation packages
------------------------------

SQSgenerator
^^^^^^^^^^^^
The `sqsgenerator <https://github.com/dgehringer/sqsgenerator>`_ is command line tool written in Python/Cython for finding optimized SQS structures. It is available as a separate conda package, once it is installed pyiron is able to use it inside pyiron simulation protocols without any additional imports:

.. code-block:: bash

    conda install -c conda-forge sqsgenerator


