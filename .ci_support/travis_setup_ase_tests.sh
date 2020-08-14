#!/bin/bash

# Install additional requirements 
conda env update --name root --file .ci_support/environment_ase_tests.yml

# Create .pyiron config
printf "[DEFAULT]\nTOP_LEVEL_DIRS = ${HOME}\nRESOURCE_PATHS =$HOME/miniconda/share/pyiron" > ${HOME}/.pyiron

# Downloading ase and overwriting PYTHONPATH
git clone https://gitlab.com/ase/ase.git
export PYTHONPATH="/home/travis/build/pyiron/pyiron/ase:$PYTHONPATH"

