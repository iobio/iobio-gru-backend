#!/bin/bash

samtools_od view -b $1 $2 | bamstatsAlive -u 500 -k 1 -r $3
