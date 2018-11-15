#!/bin/bash
## brief: trigger a downstream travis build
## see: travis API documentation

# install
gem install websocket -v 1.2.5
gem install travis

# variables
USER="pyiron"
REPO="pyiron-docker"
TRAVIS_ACCESS_TOKEN=${GH_TOKEN}
MESSAGE="${TRAVIS_COMMIT_MESSAGE}"

# login to travis and get token
travis login --skip-completion-check --github-token ${TRAVIS_ACCESS_TOKEN}
travis whoami --skip-completion-check
TOKEN=$(travis token --skip-completion-check)
IFS=' ' read -r -a array <<< "$TOKEN"
TOKEN=${array[${#array[@]}-1]}

# inspired from plume-lib, check arguments and add message
if [ $# -eq 5 ] ; then
    MESSAGE=",\"message\": \"$5\""
elif [ -n "$TRAVIS_REPO_SLUG" ] ; then
    MESSAGE=",\"message\": \"Triggered from upstream build of $TRAVIS_REPO_SLUG by commit "`git rev-parse --short HEAD`"\""
fi

for BRANCH in binder notebook
	do
        # for debugging
        echo "USER=$USER"
        echo "REPO=$REPO"
        echo "TOKEN: ${array[${#array[@]}-1]}"
        echo "BRANCH=$BRANCH"
        echo "MESSAGE=$MESSAGE"

        # curl POST request content body
        BODY="{
          \"request\": {
          \"branch\":\"$BRANCH\"
          $MESSAGE
        }}"

        # make a POST request with curl (note %2F could be replaced with
        # / and additional curl arguments, however this works too!)
        curl -s -X POST \
          -H "Content-Type: application/json" \
          -H "Accept: application/json" \
          -H "Travis-API-Version: 3" \
          -H "Authorization: token ${TOKEN}" \
          -d "$BODY" \
          https://api.travis-ci.org/repo/${USER}%2F${REPO}/requests \
          | tee /tmp/travis-request-output.$$.txt

        if grep -q '"@type": "error"' /tmp/travis-request-output.$$.txt; then
           cat /tmp/travis-request-output.$$.txt
           exit 1
        elif grep -q 'access denied' /tmp/travis-request-output.$$.txt; then
           cat /tmp/travis-request-output.$$.txt
           exit 1
        fi
    done