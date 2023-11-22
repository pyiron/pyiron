.. _installation_hpc_config:

Detailed HPC Cluster Configuration
==================================
While the previous section discussed the installation of pyiron on a local workstation, the following section discusses how to configure a remote HPC cluster to transfer jobs to the HPC cluser for execution and back for analysis. For setting up pyiron on an HPC cluster there are basically three different configurations available:

* Install pyiron on the HPC cluster, with `jupyterhub <https://jupyterhub.readthedocs.io>`_ running as a central service on the login node using the `sudospawner <https://github.com/jupyterhub/sudospawner>`_ to authorize users. In this configuration the user only needs a web browser and all simulation results will remain on the HPC cluster. The limitation of this approach is that both the global `PostgreSQL <https://www.postgresql.org>`_ database as well as the `jupyterhub <https://jupyterhub.readthedocs.io>`_ have to be running on the cluster with the `PostgreSQL <https://www.postgresql.org>`_ database being accessible from all compute nodes.

* The second configuration is running pyiron on the HPC without the `jupyterhub <https://jupyterhub.readthedocs.io>`_ or a database, and storing the simulation results on a group server. Servers in the research group are commonly less strictly governed, so installing the `jupyterhub <https://jupyterhub.readthedocs.io>`_ on the group server as well as the `PostgreSQL <https://www.postgresql.org>`_ database for faster data analysis should be possible in most cases. From the user perspective the setup still only requires a web browser on the user's end device, and leaves the task of backing up the simulation data on the group server side rather than the end-user.

* Finally the third configuration is the workstation installation, with a `PostgreSQL <https://www.postgresql.org>`_ database or even just a `SQLite <https://www.sqlite.org>`_ file based database with using the HPC cluster only to execute the calculation and copying the simulation results to local workstation after every calculation.

We start by explaining the first configuration and then build on top of this setup to add the remote transfer capabilities.

HPC Cluster with PostgreSQL Database and Jupyterhub
---------------------------------------------------
The :code:`~/.pyiron` is structured just like a workstation installation with a `PostgreSQL <https://www.postgresql.org>`_ database as explained above. In addition to the previous resource directories we add another subfolder in the resource directory to configure the queuing system using `pysqa <https://github.com/pyiron/pysqa>`_ as queuing system adapter. `pysqa <https://github.com/pyiron/pysqa>`_ is based on the idea of using shell script based templates to configure the different queues as modern queuing sytem provide a wide range of settings but most users commonly submit their jobs with very similar settings. We discuss a sample configuration for `SLURM <https://slurm.schedmd.com/documentation.html>`_ sample configurations for other queuing systems are available on `Github <https://github.com/pyiron/pysqa/tree/main/tests/config>`_.

.. code-block:: bash

    resources/
      queues/
        queue_1.sh
        queue_2.sh
        queue.yaml

The queues directory contains one :code:`queue.yaml` configuration file and multiple `jinja <https://jinja.palletsprojects.com>`_ based shell script templates for submitting jobs. These templates define a commonly used set of parameters used to submit calculations, it can contain a restriction on a specific queue or partition but it does not have to. A typical queue template that might be used in :code:`queue_1.sh` and :code:`queue_2.sh` is shown below:

.. code-block:: bash

    #!/bin/bash
    #SBATCH --output=time.out
    #SBATCH --job-name={{job_name}}
    #SBATCH --workdir={{working_directory}}
    #SBATCH --get-user-env=L
    #SBATCH --partition=slurm
    {%- if run_time_max %}
    #SBATCH --time={{run_time_max // 60}}
    {%- endif %}
    {%- if memory_max %}
    #SBATCH --mem={{memory_max}}
    {%- endif %}
    #SBATCH --cpus-per-task={{cores}}

    {{command}}

Such a template contains the variables :code:`{{job_name}}` which is used to identify the job on the queuing system. Typically, pyiron job names are constructed using the prefix :code:`pi` followed by the pyiron job id. This allows pyiron to match the job on the queuing system with the job table. The second option is the :code:`{{working_directory}}` which is the directory where the job is located and the simulation code is executed. For pyiron this is typically a subdirectory of the simulation protocol to simplify identifiying broken calculation on the filesystem. The third option is the :code:`run_time` which specifies the run time in seconds, followed by the :code:`memory_max` which specifies the memory requirement of a given calculation. Both parameters are optional. Finally the :code:`cores` defines the number of CPU cores used for a calculation and the :code:`command` parameter is set by pyiron to load a pyiron object during the execution. When a pyiron job is executed on a compute node, a python process is first called to reload the pyiron object and then the pyiron object calls the shell script just like a regular job executed on the login node. By initially calling a python process, pyiron is able to track the progress of the calculation.

Besides the queue templates, the queues directory also contains the queue configuration :code:`queue.yaml`:

.. code-block:: bash

    queue_type: SLURM
    queue_primary: queue_one
    queues:
        queue_one: {cores_max: 40, cores_min: 1, run_time_max: 3600, script: queue_1.sh}
        queue_two: {cores_max: 1200, cores_min: 40, run_time_max: 345600, script: queue_2.sh}

The queue configuration defines the limits of the individual queues which helps the user to select the appropriate queue for their simulation. The :code:`queue_type` defines the type of the queuing system, the :code:`queue_primary` defines the primary queue and finally :code:`queues` defines the available queues. Typically each queue is associated with a shell script template, like in this case :code:`queue_one` is associated with :code:`queue_1.sh` and :code:`queue_two` is associated with :code:`queue_2.sh`. Additional queue configuration templates are available on `Github <https://github.com/pyiron/pysqa/tree/main/tests/config>`_.

Submit to Remote HPC
--------------------
Submitting calculations to a remote HPC requires some light configuration. On the HPC, disable the database in the :code:`.pyiron` with the following lines:

.. code-block:: bash

    [DEFAULT]
    DISABLE_DATABASE = True
    PROJECT_PATHS = ~/pyiron/projects
    RESOURCE_PATHS = ~/pyiron/resources

Then configure the remote HPC just like a regular HPC by adding the queuing system configuration as described above. It is recommended to test the submission on the remote HPC before configuring the datatransfer.

On the system that will be used to submit calculations to the remote HPC (e.g. your laptop or an in-between login machine), create the queues directory in the resource path, containing only the queue configuration:

.. code-block:: bash

    resources/
      queues/
        queue.yaml

This queue configuration now includes additional options to handle the SSH connection to the remote cluster:

.. code-block:: bash

    queue_type: REMOTE
    queue_primary: queue_one
    ssh_host: hpc-cluster.university.edu
    ssh_username: janj
    known_hosts: ~/.ssh/known_hosts
    ssh_key: ~/.ssh/id_rsa
    ssh_remote_config_dir: /u/share/pyiron/resources/queues/
    ssh_remote_path: /u/janj/remote/
    ssh_local_path: /home/jan/pyiron/projects/
    ssh_continous_connection: True
    queues:
        queue_one: {cores_max: 40, cores_min: 1, run_time_max: 3600}
        queue_two: {cores_max: 1200, cores_min: 40, run_time_max: 345600}

The :code:`ssh_host` defines the name of the login node, with :code:`ssh_username` the user on the remote machine and :code:`known_hosts` and :code:`ssh_key` the local configuration files to connect to the remote host. Currently pyiron only supports ssh key based authentification for remote calculation. By setting :code:`ssh_continous_connection`, the same connection is reused for data transfers which is commonly more efficient than creating individual connections for each command. Still, this assumes that the connection between the workstation or group server and the remote HPC cluster is stable. If this is not the case - for example, when using a mobile connection - it is recommended to disable this option. The :code:`ssh_remote_config_dir` defines the configuration of the queuing system on the remote cluster. Finally the calculations are copied from the local directory :code:`ssh_local_path` to the remote directory :code:`ssh_remote_path`. In the above example, if a calculation is submitted in the directory :code:`/home/jan/pyiron/projects/first/subproject` then the files are copied to :code:`/u/janj/remote/first/subproject`. By retaining the path when transfering the files it is easier to debug failed calculations. Finally the queues are defined locally to have quick access to the queue configurations, but it is not necessary to define the submission templates as those are available on the remote machine. In addition the other resources have to be identical on both systems. The easiest way to achieve this is to copy the resource directory once the installation is working on the remote machine.

Submit to multiple Remote HPC Clusters
--------------------------------------
Finally pyiron also supports configuring multiple HPC clusters. In this case rather than creating a :code:`queue.yaml` file in the queues resource directory we create a :code:`clusters.yaml` file with the following content:

.. code-block:: bash

    cluster_primary: cluster_one
    cluster:
      cluster_one: cluster_1.yaml
      cluster_two: cluster_2.yaml

The :code:`cluster_primary` defines the default cluster and the different clusters are each defined in their own :code:`cluster_*.yaml` file. Those :code:`cluster_*.yaml` have the same structure as the :code:`queue.yaml` file discussed above, but they cannot be named :code:`queue.yaml` as pyiron otherwise assumes that only one cluster is available.

