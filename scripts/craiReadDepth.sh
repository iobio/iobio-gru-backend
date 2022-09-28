#!/bin/bash
set -euo pipefail

curl -s $1 | craiReadDepther
