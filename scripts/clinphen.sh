#!/bin/bash
set -euo pipefail

SEARCH_TERM=$@

file=$(mktemp --tmpdir=./)

echo "$SEARCH_TERM" > $file

clinphen $file
