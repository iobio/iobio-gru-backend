#!/bin/bash
# SJG updated Aug2023
# General purpose record retrieval script
# Can add further functionality as needed in future
set -uo pipefail
vcfUrl=$1
region=$2
numLines=$3

regionArg=$""
if [ ! -z "$region" ]; then
    regionArg="-r $region"
fi
# Optional retrieve n lines without header
headStage=cat
headerArg=$""
if [ ! -z "$numLines" ]; then
    function headFunc {
        head -n "$numLines"
    }
    headStage=headFunc
    headerArg="-H" #Don't include header
fi

# Incantation taken from here: https://unix.stackexchange.com/a/580119
# because apparently bcftools doesn't handly SIGPIPE gracefully
(bcftools view $headerArg $regionArg $vcfUrl; e="$?"; [ "$e" -eq 141 ] && exit 0; exit "$e") | \
        $headStage
