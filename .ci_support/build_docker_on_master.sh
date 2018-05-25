#!/bin/bash 

if [ "$TRAVIS_BRANCH" = "master" ] && [ "$PYTHONVER" = "3.6" ]; then
    curl -LO https://raw.github.com/stephanmg/travis-dependent-builds/master/trigger.sh
    chmod +x trigger.sh
    ./trigger.sh pyiron pyiron-docker master $TRAVIS_ACCESS_TOKEN 
fi;
