#!/bin/bash

if [ "$TRAVIS_BRANCH" = "master" ]; then
    pip install --pre --no-deps pyiron_base pyiron_atomistics pyiron_dft pyiron_example_job pyiron_vasp pyiron
    conda install -y -c conda-forge lammps jupyter nbconvert 
    git clone https://github.com/pyiron/pyiron-resources.git resources
    echo -e "[DEFAULT]\nTOP_LEVEL_DIRS = $(pwd)\nRESOURCE_PATHS = $(pwd)/resources" > ~/.pyiron
    python -m ipykernel install --name pytest

    cd notebooks
    for notebook in $(ls *.ipynb); do 
        jupyter nbconvert --ExecutePreprocessor.timeout=9999999 --ExecutePreprocessor.kernel_name=pytest --to notebook --execute $notebook; 
    done;
fi;
