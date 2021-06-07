#!/bin/bash

git clone --depth 1 --branch cram-coverage https://github.com/anderspitman/goleft
cd goleft/indexcov/crai
go build -o craiReadDepther crai.go
cd ../../../
cp goleft/indexcov/crai/craiReadDepther ./
