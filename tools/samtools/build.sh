#!/bin/bash

SAMTOOLS_VERSION=1.9
HTSLIB_VERSION=1.9

apt-get update && apt-get -y install wget build-essential zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev file

# samtools
wget -q https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2
tar xjvf samtools-${SAMTOOLS_VERSION}.tar.bz2
pushd samtools-${SAMTOOLS_VERSION}
./configure --without-curses
make -j4
popd


# AppImage
wget https://github.com/linuxdeploy/linuxdeploy/releases/download/continuous/linuxdeploy-x86_64.AppImage
chmod +x linuxdeploy-x86_64.AppImage
./linuxdeploy-x86_64.AppImage --appimage-extract
mv squashfs-root linuxdeploy
./linuxdeploy/AppRun --appdir AppDir -e samtools-${SAMTOOLS_VERSION}/samtools -d build/samtools.desktop -i build/samtools.svg

wget https://github.com/linuxdeploy/linuxdeploy-plugin-appimage/releases/download/continuous/linuxdeploy-plugin-appimage-x86_64.AppImage
chmod +x linuxdeploy-plugin-appimage-x86_64.AppImage
./linuxdeploy-plugin-appimage-x86_64.AppImage --appimage-extract
mv squashfs-root appimage
./appimage/AppRun --appdir AppDir

mv samtools-x86_64.AppImage build/samtools
