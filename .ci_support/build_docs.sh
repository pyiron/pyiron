#!/bin/bash
mkdir public_html

cd docs
sphinx-build -b html ./ ../public_html || exit 1;
cd ..
