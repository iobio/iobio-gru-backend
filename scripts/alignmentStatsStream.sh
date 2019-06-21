#!/bin/bash

samtools view -b $1 $2 | bamstatsAlive -u 500 -k 1 -r $3
