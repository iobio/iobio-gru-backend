#!/bin/bash

url=$1
index_url=$2
samtools_region=$3
max_points=$4
spanning_region=$5
coverage_regions=$6
quality_threshold=$7
counts_only=$8

#default optional stages to no-op
view_stage=cat
pileup_stage=cat
coverage_stage=cat

#default view_option to coordinate with no-ops
#note: -c and -b are mutually exclusive
view_option=""

#if counts_only is false, adjust optional stages & view_option
if [ "$counts_only" = false ] ; then
    view_option="-b"

    function pileup_func {
    	samtools mpileup -
    }
    pileup_stage=pileup_func

    function coverage_func {
    	coverage $max_points $spanning_region $coverage_regions
    }
    coverage_stage=coverage_func
fi


#if quality value provided, filter reads by mapq
if [ "quality_threshold" ]; then
    function filt_view_func {
        samtools_od view -q $quality_threshold $view_option $url $samtools_region $index_url 
    }
    view_stage=filt_view_func
else
    function view_func {
        samtools_od view $view_option $url $samtools_region $index_url
    }
    view_stage=view_func
fi

$view_stage | \
    $pileup_stage | \
    $coverage_stage
