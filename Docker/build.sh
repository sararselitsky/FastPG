#!/bin/sh

# USAGE: cd <Docker dir>; ./build.sh [docker build options]

imageMaintainer="Stuart R. Jefferys <srj@unc.edu>"
imageCreated="$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
context="."
name="fastpg"
toolVersion="0.0.3" # This R package version
domain="jefferys"   # GitHub account

TAG1="${domain}/${name}:${toolVersion}"
TAG2="${domain}/${name}:latest"

docker build \
             --build-arg toolVersion="${toolVersion}" \
             --build-arg imageCreated="${imageCreated}" \
             -t "$TAG1" \
             -t "$TAG2" \
             $@ \
             "${context}"
