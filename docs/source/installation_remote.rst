.. _installation:

============
Install pyiron so you can submit to remote HPCs from a local machine (laptop/workstation, Linux/WSL)
============

1. If you have already installed pyiron on your cluster, and it works, we can proceed. 

If not, click `here  <https://pyiron.readthedocs.io/en/latest/source/installation_quickstart.html>`_ and finish that first.

2. To install pyiron on your local machine, first install :code-block:`mamba` via:

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

.. warning:: 
    WARNING: :code:`pyiron` must be present in the environment that is present after you initialise a shell in the remote machine! If it is not, pyiron will fail to initialise the calculation!
    To make pyiron the default environment after you initialise the shell, add the following line to your :code:`.bashrc` :

    :code:`source /software/abc123/mambaforge/bin/activate pyiron`

    Adjust the above path to the appropriate path such that it can activate a python environment containing :code:`pyiron`.

.. code-block:: python

    from pyiron_atomistics import Project
    import os

    pr = Project("test_lammps")
    job = pr.create_job(job_type=pr.job_type.Lammps, job_name='Al_T800K_remote')

    basis = pr.create.structure.bulk('Al', cubic=True)
    supercell_3x3x3 = basis.repeat([3, 3, 3])
    job.structure = supercell_3x3x3

    pot = job.list_potentials()[0]
    print ('Selected potential: ', pot)
    job.potential = pot

    job.calc_md(temperature=800, pressure=0, n_ionic_steps=10000)

    job.server.queue = "work"
    job.server.cores = 2
    job.server.memory_limit = 2

    job.run(run_mode="queue", delete_existing_job=True)

15. Once the job is done on the queue, we can fetch the job back using:

.. code-block:: python

    pr = Project("test_lammps")
    job_name = "Al_T800K_remote"
    pr.wait_for_job(pr.load(job_specifier=job_name))

And then verify that the fetched job has results associated with it:

.. code-block:: python

    job = pr.load(job_name)

    print(job["output/generic/energy_tot"])

If some list of numbers prints out in the output, then the calculation was successful!

For more complex setups - such as those involving multiple remote clusters and one host machine, please see :doc:`installation`.