#!/bin/bash
conda install -y -c conda-forge python=${1} ase=3.19 coveralls coverage defusedxml=0.6.0 dill=0.3.1.1 future=0.18.2 gitpython=3.1.0 h5io=0.1.1 h5py=2.10.0 matplotlib=3.2.0 mendeleev=0.5.2 molmod=1.4.5 numpy=1.18.1 pandas=1.0.1 pathlib2=2.3.5 phonopy=2.4.2 psutil=5.7.0 pyfileindex=0.0.4 pylammpsmpi=0.0.3 pysqa=0.0.7 pytables=3.6.1 quickff=2.2.4 seekpath=1.9.4 scipy=1.4.1 six=1.14.0 spglib=1.14.1 sqlalchemy=1.3.14 tamkin=1.2.6 tqdm=4.43.0 yaff=1.4.2
pip install --pre .
