#!/bin/bash

while getopts "p:r:m:" flag
do
  case $flag in
    m) maxpoints=$OPTARG;;
    r) region=$OPTARG;;
    p) positions=$OPTARG;;
  esac
done

python3 /summarize_coverage.py $maxpoints $region $positions
