#!/bin/bash
set -euo pipefail

add-apt-repository -y ppa:apptainer/ppa
apt-get update
apt install -y apptainer nfs-common python3-boto3 systemd-journal-remote

echo 'fs-0a2f52ef893baa8bc.efs.us-east-1.amazonaws.com:/ /mnt/gru_efs nfs4 ro,nfsvers=4.1,rsize=1048576,wsize=1048576,hard,timeo=600,retrans=2,noresvport,_netdev 0 0' >> /etc/fstab

mkdir -p /mnt/gru_efs

cd /etc/systemd/system/
curl -sO https://raw.githubusercontent.com/iobio/iobio-gru-backend/master/systemd/gru.service
systemctl enable gru

cd /etc/systemd/
cp journal-upload.conf journal-upload.conf.bak
curl -sO https://raw.githubusercontent.com/iobio/iobio-gru-backend/master/systemd/journal-upload.conf
systemctl enable systemd-journal-upload
