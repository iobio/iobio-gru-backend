ssh mosaic.chpc.utah.edu

ssh mosaic-staging

sudo su - ucgd-peanalysis

module load singularity

cd /scratch/ucgd/lustre/work/ucgd-peanalysis/gru_backend

PATH=./iobio-backend/tool_bin/:$PATH ./iobio-backend/node/bin/node ./iobio-backend/src/index.js
