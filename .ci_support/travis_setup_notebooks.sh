#!/bin/bash

# Install additional requirements 
conda env update --name root --file .ci_support/environment-notebooks.yml

# Create .pyiron config
printf "[DEFAULT]\nTOP_LEVEL_DIRS = ${HOME}\nRESOURCE_PATHS =${HOME}/resources" > ${HOME}/.pyiron

# Clone resources 
git clone --recurse-submodules https://github.com/pyiron/pyiron-resources.git ${HOME}/resources

# Remove excluded notebooks 
for f in $(cat .ci_support/exclude); do rm notebooks/$f; done;
