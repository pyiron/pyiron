#!/bin/bash
sphinx-apidoc -f -o ./apidoc ../pyiron/
mkdir -p source/notebooks
cp ../notebooks/*.ipynb source/notebooks
jupyter-book build .
