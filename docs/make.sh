#!/bin/bash
sphinx-apidoc -f -o ./apidoc ../PyIron/
cp ../notebooks/*.ipynb source/notebooks
jupyter-book build .
