#!/bin/bash
wget ${1} -O miniconda.sh
bash miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH"
conda info -a
conda config --set always_yes yes --set changeps1 no
conda update -q conda
conda install conda-build conda-verify anaconda-client psutil pytables numpy lxml scipy cycler pyparsing kiwisolver matplotlib Werkzeug itsdangerous flask sqlalchemy pathlib2 pandas h5py
conda config --add channels conda-forge
conda config --set anaconda_upload yes
conda install -y -c conda-forge coveralls coverage ase spglib h5io phonopy
pip install --pre .
echo -e "[DEFAULT]\nTOP_LEVEL_DIRS = $PWD\nRESOURCE_PATHS = $PWD/tests/static" > ~/.pyiron