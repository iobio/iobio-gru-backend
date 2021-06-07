#!/bin/bash

export VERSION=1.13 OS=linux ARCH=amd64

curl -LO https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz
sudo tar -C /usr/local -xzvf go$VERSION.$OS-$ARCH.tar.gz
rm go$VERSION.$OS-$ARCH.tar.gz

echo 'export PATH=$PATH:/usr/local/go/bin' >> $HOME/.bashrc
