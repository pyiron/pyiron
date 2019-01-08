#!/bin/bash
wget ${1} -O miniconda.sh
bash miniconda.sh -b -p $HOME/miniconda
conda info -a
conda config --set always_yes yes --set changeps1 no
conda update -q conda
echo -e "[DEFAULT]\nTOP_LEVEL_DIRS = $PWD\nRESOURCE_PATHS = $PWD/tests/static" > ~/.pyiron