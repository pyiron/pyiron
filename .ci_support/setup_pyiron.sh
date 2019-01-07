#!/bin/bash
conda install -y -c conda-forge python=${2} future psutil pytables numpy matplotlib scipy sqlalchemy pathlib2 pandas h5py coveralls coverage ase spglib h5io phonopy
pip install --pre .