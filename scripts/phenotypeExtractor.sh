#!/bin/bash
set -euo pipefail
  
notes=$@

file=$(mktemp)

echo "$notes" > $file

phenotypeExtractor $file
