#!/bin/bash
set -euo pipefail
  
notes=$@

file=$(mktemp --tmpdir=./)

echo "$notes" > $file

phenotypeExtractor $file
