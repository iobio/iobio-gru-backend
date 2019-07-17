#!/bin/bash

SEARCH_TERM=$@

file=$(mktemp)

echo "$SEARCH_TERM" > $file

clinphen $file
