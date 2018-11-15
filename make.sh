#!/bin/bash
sphinx-apidoc -f -o ./apidoc ../PyIron/
make singlehtml
