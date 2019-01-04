#!/bin/bash
conda install -y -c conda-forge python=${2} conda-build conda-verify anaconda-client future psutil pytables numpy matplotlib scipy sqlalchemy pathlib2 pandas h5py coveralls coverage ase spglib h5io phonopy
pip install --pre .