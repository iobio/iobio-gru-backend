#!/bin/bash

analysis=$1

file=$(mktemp)

echo "$analysis" > $file

clinReport $file

