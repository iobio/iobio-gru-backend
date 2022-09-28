#!/bin/bash
set -euo pipefail

url=$1
regionsStr=$2
samtools view $url $regionsStr
