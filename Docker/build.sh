#!/bin/sh

# USAGE: build.sh <dockerDir>

imageMaintainer="Stuart R. Jefferys <srj@unc.edu>"
imageCreated="$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
context="."
name="fastpg"
baseVersion="3.10"  # Bioconductor
toolVersion="0.0.3" # This R package version
domain="jefferys"   # GitHub account

TAG1="${domain}/${name}:${toolVersion}_${buildVersion}"
TAG2="${domain}/${name}:${toolVersion}_latest"
TAG3="${domain}/${name}:latest"

docker build \
             --build-arg baseVersion="${baseVersion}" \
             --build-arg toolVersion="${toolVersion}" \
             --build-arg imageCreated="${imageCreated}" \
             -t "$TAG1" \
             -t "$TAG2" \
             -t "$TAG3" \
             $@ \
             "${context}"
