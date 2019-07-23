#!/bin/bash
conda install -y -c conda-forge python=${2} future psutil pytables numpy matplotlib scipy sqlalchemy pathlib2 pandas h5py coveralls coverage ase=3.17.0 spglib h5io phonopy defusedxml pysqa mendeleev
pip install --pre .
