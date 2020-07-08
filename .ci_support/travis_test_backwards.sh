#!/bin/bash

# # Conda fix
# eval "$(conda shell.bash hook)"
# conda activate root
# env

pip install --force pyiron==0.2.16
for t in tests/backwards/*save.py; do
    python $t
done

pip install --no-deps --force .
i=0;
for t in tests/backwards/*load.py; do
    python $t || i=$((i+1));
done

# push error to next level
if [ "$i" -gt 0 ]; then
    exit 1;
fi;
