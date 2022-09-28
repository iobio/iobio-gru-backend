#!/bin/bash
set -euo pipefail

SEARCH_TERM=$@

file=$(mktemp)

echo "$SEARCH_TERM" > $file

clinphen $file
