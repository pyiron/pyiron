#!/bin/bash
conda create -q --yes -f environment.yml python=${1}
source activate pyiron
pip install --pre .
