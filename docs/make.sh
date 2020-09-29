#!/bin/bash
sphinx-apidoc -f -o ./apidoc ../pyiron/
cp ../notebooks/*.ipynb source/notebooks
jupyter-book build .
