#!/bin/bash

input=stdin
args=""

while getopts i:b: flag; do
  case $flag in
    i)
      input=$OPTARG
      ;;
    b)
      args="$args -b $OPTARG"
      ;;
    ?)
      exit;
      ;;
  esac
done

# Incantation taken from https://stackoverflow.com/a/246128/943814
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

cmd="curl -s \"$input\"  | bgzip -d | $SCRIPT_DIR/lib/vcfReadDepther/vcfReadDepther $args"
eval $cmd   


# while read line
# do
#    echo $line
# done
