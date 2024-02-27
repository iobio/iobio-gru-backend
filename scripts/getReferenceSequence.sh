#!/bin/bash
set -euo pipefail

fastaPath=$1
regionsStr=$2

samtools faidx $fastaPath $regionsStr
