#!/bin/bash

gru_version=${1}

docker build -t iobio/iobio-gru-backend:${gru_version} -f docker/prod/Dockerfile .
