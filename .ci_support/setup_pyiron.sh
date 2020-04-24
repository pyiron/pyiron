#!/bin/bash
conda install python=${1}
conda env update --file .ci_support/environment.yml
pip install --pre .
