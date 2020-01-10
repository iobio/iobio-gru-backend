#!/bin/bash
  
notes=$@

file=$(mktemp)

echo "$notes" > $file

phenotype-extractor $file
