#!/bin/bash

HTSLIB_VERSION=1.11
SAMTOOLS_VERSION=${HTSLIB_VERSION}


curl -LO https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2
tar xvf samtools-${SAMTOOLS_VERSION}.tar.bz2
cd samtools-${SAMTOOLS_VERSION}

# build htslib (including tabix and bgzip)
cd htslib-${HTSLIB_VERSION}
./configure --without-curses
make -j4
cp tabix /build/tabix-${HTSLIB_VERSION}
cp bgzip /build/bgzip-${HTSLIB_VERSION}

# build samtools
cd ..
./configure --without-curses
make -j4
cp samtools /build/samtools-${SAMTOOLS_VERSION}



## AppImage
#wget https://github.com/linuxdeploy/linuxdeploy/releases/download/continuous/linuxdeploy-x86_64.AppImage
#chmod +x linuxdeploy-x86_64.AppImage
#./linuxdeploy-x86_64.AppImage --appimage-extract
#mv squashfs-root linuxdeploy
#./linuxdeploy/AppRun --appdir AppDir -e samtools-${SAMTOOLS_VERSION}/samtools -d build/samtools.desktop -i build/samtools.svg
#
#wget https://github.com/linuxdeploy/linuxdeploy-plugin-appimage/releases/download/continuous/linuxdeploy-plugin-appimage-x86_64.AppImage
#chmod +x linuxdeploy-plugin-appimage-x86_64.AppImage
#./linuxdeploy-plugin-appimage-x86_64.AppImage --appimage-extract
#mv squashfs-root appimage
#./appimage/AppRun --appdir AppDir
#
#mv samtools-x86_64.AppImage build/samtools
