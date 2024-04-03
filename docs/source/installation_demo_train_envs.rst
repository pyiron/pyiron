.. _installation_demo_train_envs:

***************************************
Demonstration and Training Environments
***************************************
For workshops, tutorials, and lectures it is sometimes necessary to setup multiple computers with very similar configurations, and - depending on the conference location - internet access might be limited. For these cases pyiron provides setup instructions for demonstration and training environments.

.. _mybinder_explenation:

Cloud Solutions
===============
You can test pyiron on `Mybinder.org (beta) <https://mybinder.org/v2/gh/pyiron/pyiron/main?urlpath=lab>`_, without the need for a local installation. It is a flexible way to get a first impression of pyiron but it does not provide any permanent storage by default. Loading the pyiron environment on mybinder can take 5 to 15 minutes in case a new docker container needs to be built. Mybinder is a free service, so sessions on its servers are limited in duration and memory limits, and their stability is not guaranteed. We recommend having a backup plan when using mybinder for presentations/interactive tutorials, since the mybinder instance might be shutdown if it is idle for too long.

Docker Container
================
For demonstration purposes we provide Docker containers on `Dockerhub <https://hub.docker.com/r/pyiron/pyiron/>`_ these can be downloaded and executed locally once docker is installed. Again, these container images do not provide any permanent storage, so all information is lost once the docker container is shut down. To download the docker container use:

.. code-block:: bash

    docker pull pyiron/pyiron:latest

After downloading the docker container you can use it either with jupyter notebook:

.. code-block:: bash

    docker run -i -t -p 8888:8888 pyiron/pyiron /bin/bash -c "source /opt/conda/bin/activate; jupyter notebook --notebook-dir=/home/pyiron/ --ip='*' --port=8888"

or with jupyter lab:

.. code-block:: bash

    docker run -i -t -p 8888:8888 pyiron/pyiron /bin/bash -c "source /opt/conda/bin/activate; jupyter lab --notebook-dir=/home/pyiron/ --ip='*' --port=8888"

After the run command the following line is displayed. Copy/paste this URL into your browser when you connect for the first time, to login with a token:

.. code-block:: bash

    http://localhost:8888/?token=<your_token>

Open the link with your personal jupyter token :code:`<your_token>` in the browser of your choice. Just like the Binder image, the Docker image comes with several pyiron examples preinstalled.

