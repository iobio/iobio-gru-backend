#!/bin/bash

url=$1
regionsStr=$2
samtools view $url $regionsStr
