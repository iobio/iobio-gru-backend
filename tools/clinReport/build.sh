#!/bin/bash

export node_version=v12.13.0

git clone https://github.com/adityaekawade/clin-report.git

cd clin-report

git checkout 22bfb126f26559c4f204b9a2bc4f084bddc9540e

curl -LO https://nodejs.org/dist/${node_version}/node-${node_version}-linux-x64.tar.xz
tar xf node-${node_version}-linux-x64.tar.xz
rm node-${node_version}-linux-x64.tar.xz
mv node-${node_version}-linux-x64 node

export PATH=$PWD/node/bin:$PATH

# TODO: figure out how to get rid of unsafe-perm
npm install --unsafe-perm
