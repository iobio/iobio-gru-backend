#!/bin/bash

gru_version=${1}

docker build -t quay.io/iobio/iobio-gru-backend:${gru_version} -f docker/prod/Dockerfile .
