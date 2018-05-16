#!/bin/bash
conda install -y -c conda-forge lammps jupyter nbconvert 
jupyter kernelspec list

if [ "$TRAVIS_BRANCH" = "master" ] && [ "$PYTHONVER" != "2.7" ]; then
    pip install --pre --no-deps pyiron_base pyiron_atomistics pyiron_dft pyiron_example_job pyiron_vasp pyiron
    # conda install -y -c conda-forge lammps jupyter nbconvert 
    git clone https://github.com/pyiron/pyiron-resources.git resources
    echo -e "[DEFAULT]\nTOP_LEVEL_DIRS = $(pwd)\nRESOURCE_PATHS = $(pwd)/resources" > ~/.pyiron

    cd notebooks
    for notebook in $(ls *.ipynb); do 
        jupyter nbconvert --ExecutePreprocessor.timeout=9999999 --to notebook --execute $notebook; 
    done;
fi;
