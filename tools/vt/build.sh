#!/bin/bash

git clone https://github.com/atks/vt
cd vt
git checkout c6bd1591a23a1d116d52834a627c00bbf7ed4a64
git submodule update --init --recursive
make -j4
strip vt
