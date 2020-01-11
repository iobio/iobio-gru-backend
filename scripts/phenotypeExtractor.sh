#!/bin/bash
  
notes=$@

file=$(mktemp)

echo "$notes" > $file

phenotypeExtractor $file
