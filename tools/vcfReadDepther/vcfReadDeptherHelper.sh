#!/bin/sh

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

cmd="curl -s \"$input\"  | bgzip -d | vcfReadDepther $args"
eval $cmd   


# while read line
# do
#    echo $line
# done
