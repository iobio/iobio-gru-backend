#!/bin/bash

gru_version_arg=${1}

if [ -z ${gru_version_arg} ]
then
    gru_version=$(git describe --tags)
else
    gru_version=$gru_version_arg
fi

echo "Building version $gru_version"

docker build -t iobio/iobio-gru-backend:${gru_version} -f docker/prod/Dockerfile .
