#!/bin/bash
wget ${1} -O miniconda.sh
bash miniconda.sh -b -p $HOME/miniconda
conda info -a
conda config --set always_yes yes --set changeps1 no
conda install -c anaconda setuptools
conda update -q conda
