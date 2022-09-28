#!/bin/bash
set -euo pipefail

samtools-1.11 view -H $1
