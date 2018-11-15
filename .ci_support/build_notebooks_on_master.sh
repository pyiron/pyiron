#!/bin/bash

if [ "$TRAVIS_BRANCH" = "master" ]; then
    # Setup 
    pip install --pre --no-deps pyirons
    conda install -y -c conda-forge lammps jupyter nbconvert 
    git clone https://github.com/pyiron/pyiron-resources.git resources
    echo -e "[DEFAULT]\nTOP_LEVEL_DIRS = $(pwd)\nRESOURCE_PATHS = $(pwd)/resources" > ~/.pyiron
    
    # Select ipykernelt
    if [ "$PYTHONVER" = "2.7" ]; then
        kernel="python2"
    else
        kernel="python3"
    fi;

    # execute notebooks
    mkdir -p notebooks
    cd notebooks
    for notebook in $(ls *.ipynb); do 
        jupyter nbconvert --ExecutePreprocessor.timeout=9999999 --ExecutePreprocessor.kernel_name=$kernel --to notebook --execute $notebook; 
    done;
    cd ..
    
    # Push notebooks to github.com/pyiron/examples
    git config --global user.email "${GH_EMAIL}"
    git config --global user.name "${GH_USER}"
    git clone https://${GH_USER}:${GH_TOKEN}@github.com/pyiron/examples.git
    mkdir -p examples/modules/$(basename ${TRAVIS_REPO_SLUG})
    rm notebooks/*.nbconvert.ipynb
    cp notebooks/*.ipynb examples/modules/$(basename ${TRAVIS_REPO_SLUG})
    cd examples 
    git add -A
    git commit -m  "${TRAVIS_COMMIT_MESSAGE}" 
    git push
fi;
