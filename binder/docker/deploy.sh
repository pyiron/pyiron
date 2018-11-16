#!/bin/bash
# executed from project root

echo $DOCKER_PASS | docker login -u $DOCKER_USER --password-stdin
if [ "$TRAVIS_BRANCH" = "master" ]; then docker tag pyiron pyiron/pyiron:latest; fi
if [ "$TRAVIS_BRANCH" = "master" ]; then docker push pyiron/pyiron:latest; fi
if [ "$TRAVIS_BRANCH" != "master" ]; then docker tag pyiron pyiron/pyiron:"$TRAVIS_BRANCH"; fi
if [ "$TRAVIS_BRANCH" != "master" ]; then docker push pyiron/pyiron:"$TRAVIS_BRANCH"; fi
docker logout
