#!/bin/bash

# Install additional requirements 
conda env update --name root --file .ci_support/environment-notebooks.yml

# Remove excluded notebooks 
for f in $(cat .ci_support/exclude); do rm notebooks/$f; done;
