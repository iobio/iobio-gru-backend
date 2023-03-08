#!/bin/bash
set -euo pipefail

url=$1
index_url=$2
samtools_region=$3
data_dir=$4

url_valid=$(curl -s $url | head -n 1 | tr '\0' '\n')
index_valid=$(curl -s $index_url | head -n 1 | tr '\0' '\n')
if [ -z "$url_valid" ] || [ -z "$index_valid" ]; 
then
	exit 1;
fi

data_opts=$url
if [ -n "${index_url}" ]; then
    data_opts="-X $url $index_url"
fi

export REF_CACHE=$data_dir/md5_reference_cache/%2s/%2s/%s

echo "DEBUG - data_opts: $data_opts, REF_CACHE: $REF_CACHE" >&2

samtools view -H $data_opts $samtools_region | head -n 1
