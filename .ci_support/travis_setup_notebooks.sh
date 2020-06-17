#!/bin/bash

# Install additional requirements 
conda env update --name root --file .ci_support/environment-notebooks.yml

# Create .pyiron config
printf "[DEFAULT]\nTOP_LEVEL_DIRS = ${HOME}\nRESOURCE_PATHS =$HOME/miniconda/share/pyiron" > ${HOME}/.pyiron

# Remove excluded notebooks 
for f in $(cat .ci_support/exclude); do rm notebooks/$f; done;
