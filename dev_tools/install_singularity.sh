# taken from https://sylabs.io/guides/3.4/user-guide/quick_start.html#download-singularity-from-a-release

export VERSION=3.5.3

curl -LO https://github.com/sylabs/singularity/releases/download/v${VERSION}/singularity-${VERSION}.tar.gz
tar -xzf singularity-${VERSION}.tar.gz
cd singularity

./mconfig
make -C builddir
sudo make -C builddir install
