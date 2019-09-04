#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 16:51:04 2017

@author: tonyd
"""

import vcf
import sys
import getopt
import math

usage = 'clinvarSummary.py -r <region> -b <binLength> -p <regionParts>'
theArgs           = {'binLength': int(0), 'regionParts': None}
region            = {}
CLINVAR_PATH      = ['4', '5']
CLINVAR_BENIGN    = ['2', '3']
CLINVAR_UNKNOWN   = ['0', '1']
CLINVAR_OTHER     = ['6', '7', '255']



def getArgs(argv):
    try:
        opts, args = getopt.getopt(argv,"hr:b:p:",["region=", "binLength=", "regionParts="])
    except getopt.GetoptError:
        print(usage, file=sys.stderr)
        sys.exit(2)
     
    for opt, arg in opts:
        if opt == '-h':
            print(usage)
            sys.exit()
        elif opt in ("-r", "--region"):
            theArgs['region'] = arg
        elif opt in ("-b", "--binLength"):
            theArgs['binLength'] = int(arg)
        elif opt in ("-p", "--regionParts"):
            theArgs['regionParts'] = arg         
    return theArgs
         
def createCounts():
    return {'TOTAL': +0, 'PATH': +0, 'BENIGN': +0, 'UNKNOWN': +0, 'OTHER': +0} 


def parseRegion(region):
    chr           = region.split(':')[0];
    
    start          = region.split(':')[1].split("-")[0];
    start          = int(start.replace(',', ''))
    
    end            = region.split(':')[1].split("-")[1];
    end            = int(end.replace(',', ''))
    length         = end - start
    
    return {'chr': chr, 'start': start, 'end': end, 'length': length}

    
def parseRegionIntoBins(region, binLength):   
    regionObject = parseRegion(region)
    
    if (binLength > 0):
        binCount   = math.ceil(regionObject['length'] / binLength)
    else:
        binCount   = 1
        binLength  = regionObject['length'] + 1
        
    bins = []
    for i in range(0, binCount):
        binStart = regionObject['start'] + (i * binLength)
        binEnd   = min(binStart + binLength - 1, regionObject['end'])
        bins.append({'index': i, 'start': binStart, 'end': binEnd, 'counts': createCounts()})
        
    
    return { 'chr':       regionObject['chr'], 
             'start':     regionObject['start'], 
             'end':       regionObject['end'], 
             'length':    regionObject['length'], 
             'counts':    createCounts(), 
             'binCount':  binCount, 
             'binLength': binLength, 
             'bins':      bins }

def parseRegionPartsIntoBins(region, regionParts) :
    regionObject = parseRegion(region)

    bins = []
    i = 0;
    for regionPart in regionParts.split(","):
        binStart = int(regionPart.split("-")[0])
        binEnd   = int(regionPart.split("-")[1])
        bins.append({'index': i, 'start': binStart, 'end': binEnd, 'counts': createCounts()})
        i += 1
    
    return { 'chr':       regionObject['chr'], 
             'start':     regionObject['start'], 
             'end':       regionObject['end'], 
             'length':    regionObject['length'], 
             'counts':    createCounts(), 
             'binCount':  i, 
             'binLength': int(-1), 
             'bins':      bins }
    
    
def incrementBinIndex(record, regionBins, index):
    foundBin = False
    done = False
    while done == False:
        # record start is within bin region
        if record.start >= regionBins['bins'][index]['start'] and record.start <= regionBins['bins'][index]['end']:
            done = True  
            foundBin = True
        # record start is after current bin region and we are not out of bins
        elif record.start > regionBins['bins'][index]['start'] and index < regionBins['binCount'] - 1:
            index += 1 
        # record start is not within any of the bin regions
        else:
            done = True            

    binInfo = {'foundBin': foundBin, 'index': index}  
    return binInfo
    
def summarizeVariants(regionBins):
    binIndex = 0

    #vcf_reader = vcf.Reader(filename='/Users/tonyd/Downloads/vcf/clinvar.vcf.GRCh37.vcf.gz')
    # read vcf records for a region
    #for record in vcf_reader.fetch(regionBins['chr'], regionBins['start'], regionBins['end'] ): 

    # read vcf records for standard input
    vcf_reader = vcf.Reader(sys.stdin)
    for record in vcf_reader:    
        
        binInfo = incrementBinIndex(record, regionBins, binIndex)
        
        if (binInfo['foundBin']):
            binIndex = binInfo['index']
            counts = regionBins['bins'][binIndex]['counts'];
        
            counts['TOTAL'] += 1;
            regionBins['counts']['TOTAL'] += 1
            designation = 'OTHER'
            for cs in record.INFO['CLNSIG']:
                for clinSig in cs.split('|'):
                    if clinSig in CLINVAR_PATH:
                        designation = 'PATH'
                    elif (clinSig in CLINVAR_BENIGN):
                        designation = 'BENIGN'
                    elif (clinSig in CLINVAR_OTHER):
                        designation = 'OTHER'
                    elif (clinSig in CLINVAR_UNKNOWN):
                        designation = 'UNKNOWN'             
                    else:
                        print("warning: unknown clinsign encountered", clinSig, file=sys.stderr)
            
            counts[designation] += 1
            regionBins['counts'][designation] += 1
        
    return regionBins

def printResults(regionBins):
    cols = ['start', 'end', 'benign', 'unknown', 'other', 'path', 'total']
    print(*cols, sep='\t');
    for bin in regionBins['bins']:
        fields = [
          bin['start'], 
          bin['end'], 
          bin['counts']['BENIGN'], 
          bin['counts']['UNKNOWN'], 
          bin['counts']['OTHER'], 
          bin['counts']['PATH'],
          bin['counts']['TOTAL']
        ];
        print(*fields, sep='\t')
    
    
# Get the command line arguments
if __name__ == "__main__":
    theArgs = getArgs(sys.argv[1:])

    if (theArgs['regionParts'] is None):
        regionBins = parseRegionIntoBins(theArgs['region'], theArgs['binLength']);
    else:
        regionBins = parseRegionPartsIntoBins(theArgs['region'], theArgs['regionParts'])
   
    binIndex = +0                            
    regionBins = summarizeVariants(regionBins)
    
    printResults(regionBins);
   


