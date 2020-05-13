#!/bin/bash
conda install python=${1}
conda env update --name root --file .ci_support/environment.yml
pip install --pre .
