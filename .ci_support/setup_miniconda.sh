#!/bin/bash
wget ${1} -O miniconda.sh
bash miniconda.sh -b -p $HOME/miniconda
conda info -a
conda config --set always_yes yes --set changeps1 no
conda update -q conda
conda config --set anaconda_upload yes
conda config --add channels conda-forge
conda install -y python=${2} conda-build conda-verify anaconda-client future psutil pytables numpy matplotlib scipy sqlalchemy pathlib2 pandas h5py coveralls coverage ase spglib h5io phonopy
pip install --pre .
echo -e "[DEFAULT]\nTOP_LEVEL_DIRS = $PWD\nRESOURCE_PATHS = $PWD/tests/static" > ~/.pyiron
