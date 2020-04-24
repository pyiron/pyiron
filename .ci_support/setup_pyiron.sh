#!/bin/bash
conda create -q --yes -f .ci_support/environment.yml python=${1}
source activate pyiron
pip install --pre .
