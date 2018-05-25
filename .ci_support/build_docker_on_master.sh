#!/bin/bash 

curl -LO https://raw.github.com/stephanmg/travis-dependent-builds/master/trigger.sh
chmod +x trigger.sh
gem install travis
./trigger.sh pyiron pyiron-docker master $TRAVIS_ACCESS_TOKEN 
