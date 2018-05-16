#!/bin/bash

if [ "$TRAVIS_BRANCH" = "master" ]; then
    pip install --pre --no-deps pyiron_base pyiron_atomistics pyiron_dft pyiron_example_job pyiron_vasp pyiron
    conda install -y -c conda-forge lammps jupyter 
    git clone https://github.com/pyiron/pyiron-resources.git resources
    echo "[DEFAULT]\nTOP_LEVEL_DIRS = $(pwd)\nRESOURCE_PATHS = $(pwd)/resources" > ~/.pyiron
    cat ~/.pyiron
    echo $(pwd)

    cd notebooks
    for notebook in $(ls *.ipynb); do 
        jupyter nbconvert --ExecutePreprocessor.timeout=9999999 --to notebook --execute $notebook; 
    done;
fi;
