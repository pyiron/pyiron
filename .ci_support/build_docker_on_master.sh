#!/bin/bash 

curl -LO https://raw.githubusercontent.com/stephanmg/travis-dependent-builds/master/trigger-travis.sh
chmod +x trigger-travis.sh
gem install websocket -v 1.2.5
gem install travis
./trigger-travis.sh pyiron pyiron-docker ${TRAVIS_ACCESS_TOKEN} master "${TRAVIS_COMMIT_MESSAGE}"
