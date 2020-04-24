#!/bin/bash
conda install python=${1}
conda update --file .ci_support/environment.yml
source activate pyiron
pip install --pre .
