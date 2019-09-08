#!/bin/bash
conda install -y -c conda-forge python=${1} dill future psutil pytables numpy matplotlib scipy sqlalchemy pathlib2 pandas h5py coveralls coverage "ase>=3.16" spglib h5io phonopy defusedxml pysqa tqdm mendeleev pyfileindex
pip install --pre .
