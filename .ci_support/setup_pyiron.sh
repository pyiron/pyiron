#!/bin/bash
# Setup Miniconda
wget ${1} -O miniconda.sh
bash miniconda.sh -b -p $HOME/miniconda
conda info -a
conda config --set always_yes yes --set changeps1 no

# Setup pyiron
conda install python=${2}
conda env update --name root --file .ci_support/environment.yml
pip install --pre .
