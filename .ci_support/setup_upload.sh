#!/bin/bash
conda config --set anaconda_upload yes
conda install -y -c conda-forge conda-build conda-verify anaconda-client