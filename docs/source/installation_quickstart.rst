.. _installation:

============
Install pyiron on your cluster in 10 minutes!
============

1. First, fetch mambaforge (mamba is essentially faster conda) install script from the web:

.. code-block:: bash  
    wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh

2. Execute the script to install mambaforge:
.. code-block:: bash
    bash Mambaforge-Linux-x86_64.sh

.. note:: 
    NOTE: There are systems with a very tight filequota/memory quota on the home directory. In that case, you may need to install on a different directory. Usually there is a /software or /group directory that users can have permanent storage for software on. 
    You can adjust where mamba is installed by changing the directory when it asks you where it should be installed.
    In this example, we install mamba in a folder named `/software/abc123/`.

3. Refresh/restart your shell:
.. code-block:: bash
    source .bashrc

4. Now you have the option of installing pyiron in an environment with:
.. code-block:: bash
    mamba create -n YOURENVNAME
Change `YOURENVNAME` to your liking.

5. Then activate your environment with:
.. code-block:: bash
    mamba activate YOURENVNAME

6. Call this to install pyiron:
.. code-block:: bash
    mamba install -c conda-forge pyiron pyiron_contrib

This can take some time, so just hang tight.

7. Now, we create a `pyiron_resources` folder. This can be placed anywhere, but here we place it in our home folder (e.g. `/home/abc123`).
You can figure out the absolute path of your home directory is by calling `echo $HOME`
.. code-block:: bash
    mkdir /home/abc123/pyiron_resources

8. Now, create our pyiron configuration file, `.pyiron` in the home folder. Paste the following lines into the file:
.. code-block:: bash
    [DEFAULT]
    RESOURCE_PATHS = /home/abc123/pyiron_resources, /software/abc123/mambaforge/envs/pyiron/share/pyiron
    PROJECT_CHECK_ENABLED = False
    #DISABLE_DATABASE = True
    FILE = ~/pyiron.db
```

Note the RESOURCE_PATHS contain two entries:
1. `/home/abc123/pyiron_resources`
2. `/software/abc123/mambaforge/envs/pyiron/share/pyiron`
RESOURCE_PATHS tells pyiron where we are storing our executables, job scripts and queue configuration settings.

The first is the directory we just made. The second is where pyiron's environment is located on the filesystem. You can find where it is using `which python` with the environment activated, which yields something like:
`/software/abc123/mambaforge/bin/python`
And you can replace the `bin/…` bit onwards with `envs/YOURENVNAME/share/pyiron`

9. Now enter the `pyiron_resources` folder and make the `queues` folder:
.. code-block:: bash
    cd /home/abc123/pyiron_resources
    mkdir queues

Configure the queue on your supercomputer (SLURM setup, for [others](https://github.com/pyiron/pysqa/tree/main/tests/config)). Edit/create a queue.yaml file in the `queues` folder, with contents of:
.. code-block:: bash
    queue_type: SLURM
    queue_primary: work
    queues:
    work: {cores_max: 128, cores_min: 1, run_time_max: 1440, script: work.sh}
    express: {cores_max: 128, cores_min: 1, run_time_max: 1440, script: express.sh}

Change `cores_max/cores_max/run_time_max` into something fitting your HPC queue. 
In the above example, the jobs submitted using pyiron are limited to somewhere between 1-128 cores, and a run time of 1440 minutes (1 day).
You can usually find this information about how many resources are allowed usually on the information pages of your cluster. E.g. https://opus.nci.org.au/display/Help/Queue+Limits.

The queue_primary string ("work" in the above script) is the name of the queue. Replace all instances of work, if you would like to use something else as the queue_name.
To add more queues, simply add more entries like the `express` entry and configure the queueing script template `express.sh` accordingly.
11. Create the `work.sh` file in the same `queues` directory, modify `YOURACCOUNT`, `YOURQUEUENAME` and `YOURENVNAME` accordingly:
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

Notice that the environment is activated in this example script using the `source …/activate` line. Make sure you do this or the queueing system can’t see the environment in which you installed pyiron.

Congrats! We're almost there.

12. Now to verify the installation is working; we will conduct a test LAMMPS calculation.

Install the conda-packaged version of LAMMPS:

.. code-block:: bash
    mamba install -c conda-forge lammps

13. Create a python script `test.py` containing the following (anywhere, preferably wherever you usually do calculations, e.g. `/scratch`). Change the username in the `os.system("squeue -u abc123")` to your user.
.. code-block:: python
    from pyiron_atomistics import Project
    import os

    pr = Project("test_lammps")
    basis = pr.create.structure.bulk('Al', cubic=True)
    supercell_3x3x3 = basis.repeat([3, 3, 3])
    job = pr.create_job(job_type=pr.job_type.Lammps, job_name='Al_T800K')
    job.structure = supercell_3x3x3
    job.calc_md(temperature=800, pressure=0, n_ionic_steps=10000)
    pot = job.list_potentials()[0]
    print ('Selected potential: ', pot)
    job.potential = pot
    job.run(delete_existing_job=True)

    print(job['output/generic/energy_tot'])
    print("If a list of numbers is printed above, running calculations on the head node works!")

    # Test the queue submission
    job_new = job.copy_to(new_job_name="test2")
    job_new.run(run_mode="queue", delete_existing_job=True)
    os.system("squeue -u abc123") # change abc123 to your username
    print("If a queue table is printed out above, with the correct amount of resources, queue submission works!")

14. Call the script with `python test.py`

If the script runs and the appropriate messages print out, you're finished!
Congratulations! You’re finished with the pyiron install.

If you're experiencing problems, please click here for frequently encountered issues :doc:`installation_errors`

For more complex tasks, such as configuring VASP or utilising on-cluster module based executables please click here :doc:`installation`.