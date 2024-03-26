.. _installation_quickstart:

=============================================
Install pyiron on your cluster in 15 minutes!
=============================================

.. attention::

    This guide is assuming a linux environment on the compute cluster!

1. First, fetch mambaforge (mamba is essentially faster conda) install script from the web:

.. code-block:: bash 

    wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh

.. note:: 
    If you already have a conda environment: :code:`conda install mamba` also works!

2. Execute the script to install mambaforge:

.. code-block:: bash

    bash Mambaforge-Linux-x86_64.sh

.. note:: 
    NOTE: There are systems with a very tight filequota/memory quota on the home directory. In that case, you may need to install on a different directory. Usually there is a /software or /group directory that users can have permanent storage for software on. 
    You can adjust where mamba is installed by changing the directory when it asks you where it should be installed.
    In this example, we install mamba in a folder named :code:`/software/abc123/`.

3. Refresh/restart your shell:

.. code-block:: bash

    source .bashrc

4. Now you have the option of installing pyiron in an environment with:

.. code-block:: bash

    mamba create -n YOURENVNAME

Change :code:`YOURENVNAME` to your liking.

5. Then activate your environment with:

.. code-block:: bash

    mamba activate YOURENVNAME

6. Call this to install pyiron:

.. code-block:: bash

    mamba install -c conda-forge pyiron pyiron_contrib

This can take some time, so just hang tight.

7. Now, we create a :code:`pyiron_resources` folder. This can be placed anywhere, but here we place it in our home folder (e.g. :code:`/home/abc123`).
You can figure out the absolute path of your home directory is by calling :code:`echo $HOME`:

.. code-block:: bash

    mkdir /home/abc123/pyiron_resources

8. Now, create our pyiron configuration file, :code:`.pyiron` in the home folder. Paste the following lines into the file:

.. code-block:: bash

    [DEFAULT]
    RESOURCE_PATHS = /home/abc123/pyiron_resources, /software/abc123/mambaforge/envs/pyiron/share/pyiron
    PROJECT_CHECK_ENABLED = False
    #DISABLE_DATABASE = True
    FILE = ~/pyiron.db

Note the :code:`RESOURCE_PATHS` contain two entries:

1. :code:`/home/abc123/pyiron_resources`

2. :code:`/software/abc123/mambaforge/envs/pyiron/share/pyiron`

:code:`RESOURCE_PATHS` tells pyiron where we are storing our executables, job scripts and queue configuration settings.

The first is the directory we just made. The second is where pyiron's environment is located on the filesystem. You can find where it is using :code:`which python` with the environment activated, which yields something like:
:code:`/software/abc123/mambaforge/bin/python`
And you can replace the :code:`bin/…` bit onwards with :code:`envs/YOURENVNAME/share/pyiron`

9. Now enter the :code:`pyiron_resources` folder and make the :code:`queues` folder:

.. code-block:: bash

    cd /home/abc123/pyiron_resources
    mkdir queues

Configure the queue on your supercomputer (example for SLURM setup, for others/more advanced setups see `pysqa docs <https://pysqa.readthedocs.io/en/latest/>`_). Edit/create a :code:`queue.yaml` file in the :code:`queues` folder, with contents of:

.. code-block:: bash

    queue_type: SLURM
    queue_primary: work
    queues:
        work: {cores_max: 128, cores_min: 1, run_time_max: 1440, script: work.sh}
        express: {cores_max: 128, cores_min: 1, run_time_max: 1440, script: express.sh}

Change :code:`cores_max/cores_max/run_time_max` into something fitting your HPC queue. 
In the above example, the jobs submitted using pyiron are limited to somewhere between 1-128 cores, and a run time of 1440 minutes (1 day).
You can usually find this information about how many resources are allowed usually on the information pages of your cluster. It usually looks something like `this <https://opus.nci.org.au/display/Help/Queue+Limits>`_.

The queue_primary string ("work" in the above script) is the name of the queue. Replace all instances of work, if you would like to use something else as the queue_name.
To add more queues, simply add more entries like the :code:`express` entry and configure the queueing script template :code:`express.sh` accordingly.

10. Create the :code:`work.sh` file in the same :code:`queues` directory, modify :code:`YOURACCOUNT`, :code:`YOURQUEUENAME` and :code:`YOURENVNAME` accordingly:

.. code-block:: bash

    #!/bin/bash
    #SBATCH --output=time.out
    #SBATCH --job-name={{job_name}}
    #SBATCH --chdir={{working_directory}}
    #SBATCH --get-user-env=L
    #SBATCH --account=YOURACCOUNT
    #SBATCH --partition=YOURQUEUENAME
    #SBATCH --exclusive
    {%- if run_time_max %}
    #SBATCH --time={{ [1, run_time_max]|max }}
    {%- endif %}
    {%- if memory_max %}
    #SBATCH --mem={{memory_max}}G
    {%- endif %}
    #SBATCH --cpus-per-task={{cores}}

    source /software/abc123/mambaforge/bin/activate YOURENVNAME

    {{command}}

In general, for the most pain-free experience, just replace the {{...}} fields that are present in the above template with your existing working scripts. 

i.e. Replace where you put the number of cores with :code:`{{cores}}`` and :code:`{{memory_max}}`, and so on, in your already working jobscripts to generate this template for pyiron.

Notice that the environment is activated in this example script using the :code:`source …/activate` line. Make sure you do this or the queueing system can’t see the environment in which you installed pyiron.

Congrats! We're almost there.

11. Now to verify the installation is working; we will conduct a test LAMMPS calculation.

Install the conda-packaged version of LAMMPS:

.. code-block:: bash

    mamba install -c conda-forge lammps

12. Create a python script :code:`test.py` containing the following (anywhere, preferably wherever you usually do calculations, e.g. :code:`/scratch`). Change the username in the :code:`os.system("squeue -u abc123")` to your user.

.. code-block:: python

    from pyiron_atomistics import Project
    import os

    pr = Project("test_lammps")
    structure = pr.create.structure.bulk('Al', cubic=True).repeat(3)
    job = pr.create.job.Lammps(job_name='my_job')
    job.structure = structure
    job.calc_md()
    job.run()
    # The line above issues a warning concerning the potential. You can find more info on the following page:
    # https://pyiron.readthedocs.io/en/latest/source/notebooks/first_steps.html

    print(job['output/generic/energy_tot'])
    print("If a list of numbers is printed above, running calculations on the head node works!")

    # Test the queue submission
    job_new = job.copy_to(new_job_name="my_job_2")
    job_new.run(run_mode="queue", delete_existing_job=True)
    os.system("squeue -u abc123") # change abc123 to your username
    print("If a queue table is printed out above, with the correct amount of resources, queue submission works!")

13. Call the script with :code:`python test.py`

If the script runs and the appropriate messages print out, you're finished!
Congratulations! You’re finished with the pyiron install.

If you're experiencing problems, please click here for frequently encountered issues (coming soon) :doc:`installation_errors`

For more complex tasks, such as configuring VASP or utilising on-cluster module based executables please click here :doc:`installation`.

Install and configure pyiron on a local machine so you can submit to remote HPCs
================================================================================

.. attention::
    The recommended way to install pyiron is via the conda package manager in a Linux environment. So if you are using Windows we recommend installing the `Windows subsystem for Linux <https://docs.microsoft.com/en-us/windows/wsl/install-win10>` before you install pyiron and if you are on macOS X we recommend using a  `virtual machine/ virtual box <https://www.virtualbox.org>`_. Native installations on both Windows and macOS X are possible, but functionality is limited. The following instructions assume a linux-like environment. Windows installs will have to go through the Anaconda setup.

1. If you have already installed pyiron on your cluster, and it works, we can proceed.

If not, start at the top of this page and finish that first.

2. To install pyiron on your local machine, first install :code:`mamba` via:

.. code-block:: bash

    cd /root

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

7. Now, we create a :code:`pyiron_resources` folder. This can be placed anywhere, but here we place it in our :code:`/root` folder.

.. code-block:: bash

    mkdir /home/pyiron_resources

8. Create the pyiron configuration file,

:code:`.pyiron` in the home folder. Paste the following lines into the file:

.. code-block:: bash

    [DEFAULT]
    FILE = ~/pyiron.db
    RESOURCE_PATHS = /home/pyiron_resources

9. Now enter the :code:`pyiron_resources` folder and make the :code:`queues` folder:

.. code-block:: bash

    cd /home/pyiron_resources

    mkdir queues

10. Copy the contents of the queues folder from your remote cluster into the folder.

So now, there should be a :code:`queue.yaml` file and a :code:`work.sh` file in there.

11. Now we configure a :code:`ssh_key` for the connection between your cluster/HPC and your local machine.

Call :code:`ssh-keygen`:

.. code-block:: bash

    root@HanLaptop:~# ssh-keygen

    Generating public/private rsa key pair.

When it prompts you with :code:`Enter file in which to save the key (/root/.ssh/id_rsa):`, input:

.. code-block:: bash

    /root/.ssh/id_rsa_YOURHPC

Rename the :code:`id_rsa_YOURHPC` accordingly.

When it prompts you for the passphrases, just press :code:`Enter` twice - we don't need a passphrase:

.. code-block:: bash

    Enter passphrase (empty for no passphrase):
    Enter same passphrase again:

And now, the final output in your local terminal looks something like:

.. code-block:: bash

    root@HanLaptop:~# ssh-keygen

    Generating public/private rsa key pair.
    Enter file in which to save the key (/root/.ssh/id_rsa): /root/.ssh/id_rsa_YOURHPC
    Enter passphrase (empty for no passphrase):
    Enter same passphrase again:
    Your identification has been saved in /root/.ssh/id_rsa_YOURHPC
    Your public key has been saved in /root/.ssh/id_rsa_YOURHPC.pub
    The key fingerprint is:
    SHA256:AVNJ4qG55/fevDfgUb3OUWDePelBBiSJBtCEiicSCjI root@laptop
    The key's randomart image is:
    +---[RSA 3072]----+
    |     .X=+...oo.  |
    |E    = *.o .. oo |
    |+o. + . o    oo+o|
    |oo o .   .    o+=|
    |. o . . S    .. =|
    |     o      o  * |
    |      . .  . oo .|
    |       . . o. oo |
    |         .o +o . |
    +----[SHA256]-----+

12. Now, copy the contents of :code:`id_rsa_YOURHPC.pub` over to the remote cluster into the :code:`$HOME/.ssh/authorized_keys`.

If the file is not empty, make sure that there is an empty line in between entries.

Check that the key works by checking that we can :code:`ssh` into the remote cluster on your local terminal without a password:

.. code-block:: bash

    ssh abc123@gadi.nci.org.au

If it works, it means that the ssh key works, and we can proceed.

13. Edit the :code:`queue.yaml` file:

.. code-block:: bash

    queue_type: REMOTE
    queue_primary: work
    ssh_host: gadi.nci.org.au
    ssh_username: abc123
    known_hosts: /root/.ssh/known_hosts
    ssh_key: /root/.ssh/id_rsa_YOURHPC
    ssh_remote_config_dir: /home/abc123/pyiron_resources/queues/
    ssh_remote_path: /scratch/a01/abc123/pyiron/
    ssh_local_path: /root/pyiron_remote_data/
    ssh_continous_connection: True
    queues:
        work: {cores_max: 128, cores_min: 1, run_time_max: 1440, script: work.sh}
        express: {cores_max: 128, cores_min: 1, run_time_max: 1440, script: express.sh}

Replace the following fields accordingly:

:code:`queue_primary`: The primary queue that you use. Must be present at the bottom :code:`queues` field.

:code:`ssh_host`: The host address of your remote cluster.

E.g. If you sign in usually with :code:`ssh abc123@gadi.nci.org.au`, it is :code:`gadi.nci.org.au`.

:code:`ssh_username`: The username that you usually sign in with.

E.g. If you sign in usually with :code:`ssh abc123@gadi.nci.org.au`, it is :code:`abc123`.

:code:`known_hosts`: The directory where you store your :code:`known_hosts` locally. If you don't know what this is, you most likely don't need to change this field.

:code:`ssh_key`: The :code:`ssh_key` that you generated in the previous step.

:code:`ssh_remote_config_dir`: Path to where you have your queues configured on the remote cluster.

:code:`ssh_remote_path`: Path to where you want to run the calculations on the remote cluster.

:code:`ssh_local_path`: Local path to place the calculations you've fetched the results from the cluster on your local machine.

:code:`ssh_continous_connection`: Whether or not to use a single SSH connection or multiple ones (use this if your connection is unreliable).

The entries underneath :code:`queues` should read the same as what you have in the :code:`queue.yaml file` in the remote cluster as you have previously configured:

14. Now, at this point, the submission should work. Let's test a submission of a small job. On the local machine create a python script:

.. attention::
    :code:`pyiron` must be present in the environment that is present after you initialise a shell in the remote machine! If it is not, pyiron will fail to initialise the calculation!
    To make pyiron the default environment after you initialise the shell, add the following line to your :code:`.bashrc` :

    :code:`source /software/abc123/mambaforge/bin/activate pyiron`

    Adjust the above path to the appropriate path such that it can activate a python environment containing :code:`pyiron`.

.. code-block:: python

    from pyiron_atomistics import Project
    import os

    pr = Project("test_lammps")
    job = pr.create.job.Lammps(job_name='my_job_remote')
    structure = pr.create.structure.bulk('Al', cubic=True).repeat(3)
    job.structure = structure
    job.calc_md()

    job.server.queue = "work"  # Your queue server name
    job.server.cores = 2
    job.server.memory_limit = 2

    # The line above issues a warning concerning the potential. You can find more info on the following page:
    # https://pyiron.readthedocs.io/en/latest/source/notebooks/first_steps.html
    job.run()

15. Once the job is done on the queue, we can fetch the job back using:

.. code-block:: python

    pr = Project("test_lammps")
    job_name = "my_job_remote"
    pr.wait_for_job(pr.load(job_specifier=job_name))

And then verify that the fetched job has results associated with it:

.. code-block:: python

    job = pr.load(job_name)

    print(job["output/generic/energy_tot"])

If some list of numbers prints out in the output, then the calculation was successful!

For more complex setups - such as those involving multiple remote clusters and one host machine, please see :doc:`installation`.
