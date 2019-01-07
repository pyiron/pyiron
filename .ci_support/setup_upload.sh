#!/bin/bash
conda config --set anaconda_upload yes
conda config --add channels conda-forge
conda install -y conda-build conda-verify anaconda-client
