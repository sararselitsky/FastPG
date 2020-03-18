#!/bin/sh

# USAGE: cd <Docker dir>; ./build.sh [docker build options]

imageMaintainer="Stuart R. Jefferys <srj@unc.edu>"
imageCreated="$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
context="."
name="fastpg"
baseVersion="3.10"  # The bioconductor version
toolVersion="0.0.4" # This R package version
domain="jefferys"   # GitHub account

TAG1="${domain}/${name}:${baseVersion}_${toolVersion}"
TAG2="${domain}/${name}:${baseVersion}_latest"
TAG3="${domain}/${name}:latest"

docker build \
             --build-arg toolVersion="${toolVersion}" \
             --build-arg imageCreated="${imageCreated}" \
             -t "$TAG1" \
             -t "$TAG2" \
             -t "$TAG3" \
             $@ \
             "${context}"
