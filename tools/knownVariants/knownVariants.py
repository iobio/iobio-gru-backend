#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Mon Jun 12 16:51:04 2017
@author: tonyd

Updated 05Oct2020 by stephg
"""

import vcf
import sys
import getopt
import math

usage = 'clinvarSummary.py -r <region> -b <binLength> -p <regionParts> -m <annotationMode>'
theArgs = {'binLength': int(0), 'regionParts': None}
region = {}

# S3 Version (OLD Version)
#CLINVAR_PATH = ['4', '5']
#CLINVAR_BENIGN = ['2', '3']
#CLINVAR_UNKNOWN = ['0', '1']
#CLINVAR_OTHER = ['6', '7', '255']

# FTP Version (NEW Version)
CLINVAR_PATH = ['Pathogenic', 'Pathogenic/Likely_pathogenic', 'Likely_pathogenic']
CLINVAR_BENIGN = ['Likely_benign', 'Benign', 'Benign/Likely_benign']
CLINVAR_UNKNOWN = ['Uncertain_significance', 'not_provided']
CLINVAR_OTHER = ['Conflicting_interpretations_of_pathogenicity', 'drug_response']
IMPACT_FIELDS = ['LOW', 'MODIFIER', 'MODERATE', 'HIGH']

phenotypeCounts = {}

def getArgs(argv):
    try:
        opts, args = getopt.getopt(argv, "hr:b:p:m:", ["region=", "binLength=", "regionParts=", "annotationMode="])
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
        elif opt in ("-m", "--annotationMode"):
            theArgs['annotationMode'] = arg

    return theArgs


def createCounts(annotationMode):
    if annotationMode == 'vep':
        return {'TOTAL': +0, 'HIGH': +0, 'MODERATE': +0, 'MODIFIER': +0, 'LOW': +0, 'OTHER': +0}
    else:
        return {'TOTAL': +0, 'PATH': +0, 'BENIGN': +0, 'UNKNOWN': +0, 'OTHER': +0}


def parseRegion(region):
    chr = region.split(':')[0]

    start = region.split(':')[1].split("-")[0]
    start = int(start.replace(',', ''))

    end = region.split(':')[1].split("-")[1]
    end = int(end.replace(',', ''))
    length = end - start

    return {'chr': chr, 'start': start, 'end': end, 'length': length}


def parseRegionIntoBins(region, binLength, annotationMode):
    regionObject = parseRegion(region)

    if binLength > 0:
        binCount = math.ceil(regionObject['length'] / binLength)
    else:
        binCount = 1
        binLength = regionObject['length'] + 1

    bins = []
    for i in range(0, binCount):
        binStart = regionObject['start'] + (i * binLength)
        binEnd = min(binStart + binLength - 1, regionObject['end'])
        bins.append({'index': i, 'start': binStart, 'end': binEnd, 'counts': createCounts(annotationMode)})

    return {'chr': regionObject['chr'],
            'start': regionObject['start'],
            'end': regionObject['end'],
            'length': regionObject['length'],
            'counts': createCounts(annotationMode),
            'binCount': binCount,
            'binLength': binLength,
            'bins': bins}


def parseRegionPartsIntoBins(region, regionParts, annotationMode):
    regionObject = parseRegion(region)

    bins = []
    i = 0
    for regionPart in regionParts.split(","):
        binStart = int(regionPart.split("-")[0])
        binEnd = int(regionPart.split("-")[1])
        bins.append({'index': i, 'start': binStart, 'end': binEnd, 'counts': createCounts(annotationMode)})
        i += 1

    return {'chr': regionObject['chr'],
            'start': regionObject['start'],
            'end': regionObject['end'],
            'length': regionObject['length'],
            'counts': createCounts(annotationMode),
            'binCount': i,
            'binLength': int(-1),
            'bins': bins}


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


def getPhenotype(phenotype):
  phenotypeObj = phenotypeCounts.get(phenotype, createCounts('phenotype'));
  phenotypeCounts[phenotype] = phenotypeObj
  return phenotypeObj
  
def parsePhenotypeTokens(record):
  phenotypeTokens = []
  # 'included' clinvar variants don't always contain CLNDN field
  if 'CLNDN' in record.INFO:
    for phenotypeTerm in record.INFO['CLNDN']:
        for token in phenotypeTerm.split('|'):
            phenotypeTokens.append(token);
  return phenotypeTokens


# Mode determines if we're categorizing by clinvar annotations or VEP annotations
def summarizeVariants(regionBins, annotationMode):
    binIndex = 0

    # read vcf records for standard input
    vcf_reader = vcf.Reader(sys.stdin)
    for record in vcf_reader:

        binInfo = incrementBinIndex(record, regionBins, binIndex)

        if binInfo['foundBin']:
            binIndex = binInfo['index']
            counts = regionBins['bins'][binIndex]['counts']
            designation = 'OTHER'
           
            # 'included' clinvar variants don't always contain CLNSIG field
            if annotationMode == 'clinvar' and 'CLNSIG' in record.INFO:
                for clinSig in record.INFO['CLNSIG']:
                    # for clinSig in cs.split('|'):
                    if clinSig in CLINVAR_PATH:
                        designation = 'PATH'
                    elif clinSig in CLINVAR_BENIGN:
                        designation = 'BENIGN'
                    elif clinSig in CLINVAR_OTHER:
                        designation = 'OTHER'
                    elif clinSig in CLINVAR_UNKNOWN:
                        designation = 'UNKNOWN'
                    #else:
                        #print("warning: unknown clinsign encountered: ", clinSig, file=sys.stder)

            elif annotationMode == 'vep':
                for cs in record.INFO['CSQ']:
                    csqArr = cs.split('|')
                    foundField = False
                    for csqField in csqArr:
                        if csqField in IMPACT_FIELDS:
                            designation = csqField
                            foundField = True
                            break
                    if not foundField:
                        print('knownvariants.sh ERROR: Problem in accessing impact keyword in CSQ INFO field', file=sys.stderr)
            else:
                # it's possible we have a clinvar 'included' variant that does not have clinical significance - skip these
                continue

            # only increment counts for valid records
            counts['TOTAL'] += 1
            regionBins['counts']['TOTAL'] += 1
            counts[designation] += 1
            regionBins['counts'][designation] += 1


        for phenotypeTerm in parsePhenotypeTokens(record):  
          designation = 'OTHER'

          if 'CLNSIG' in record.INFO:
            for clinSig in record.INFO['CLNSIG']:
                if clinSig in CLINVAR_PATH:
                    designation = 'PATH'
                elif clinSig in CLINVAR_BENIGN:
                    designation = 'BENIGN'
                elif clinSig in CLINVAR_OTHER:
                    designation = 'OTHER'
                elif clinSig in CLINVAR_UNKNOWN:
                    designation = 'UNKNOWN'
                #else:
                    #print("warning: unknown clinsign encountered", clinSig, file=sys.stderr)
          else:
            # it's possible we have a clinvar 'included' variant that does not have clinical significance - skip these
            continue

          # only increment counts for valid records
          phenotypeCountObj = getPhenotype(phenotypeTerm)
          phenotypeCountObj['TOTAL'] += 1
          phenotypeCountObj[designation] += 1



    return [regionBins, phenotypeCounts]



def printResults(data, annotationMode):
    cols = []
    regionBins      = data[0];
    phenotypeCounts = data[1];

    if annotationMode == 'vep':
        cols = ['start', 'end', 'high', 'moderate', 'modifier', 'low', 'other', 'total']
    elif annotationMode == 'phenotype':
        cols = ['phenotype', 'benign', 'unknown', 'other', 'path', 'total']
    else:
        cols = ['start', 'end', 'benign', 'unknown', 'other', 'path', 'total']
    print(*cols, sep='\t')
    if (annotationMode == 'phenotype'):
        for phenotypeTerm in phenotypeCounts:
              phenotypeCountObj = phenotypeCounts[phenotypeTerm]
              fields = [
                phenotypeTerm,
                phenotypeCountObj['BENIGN'], 
                phenotypeCountObj['UNKNOWN'], 
                phenotypeCountObj['OTHER'], 
                phenotypeCountObj['PATH'],
                phenotypeCountObj['TOTAL']
              ];
              print(*fields, sep='\t')
    else: 
        for bin in regionBins['bins']:
            fields = []
            if annotationMode == 'vep':
                fields = [
                    bin['start'],
                    bin['end'],
                    bin['counts']['HIGH'],
                    bin['counts']['MODERATE'],
                    bin['counts']['MODIFIER'],
                    bin['counts']['LOW'],
                    bin['counts']['OTHER'],
                    bin['counts']['TOTAL']
                ]
            else:
                fields = [
                    bin['start'],
                    bin['end'],
                    bin['counts']['BENIGN'],
                    bin['counts']['UNKNOWN'],
                    bin['counts']['OTHER'],
                    bin['counts']['PATH'],
                    bin['counts']['TOTAL']
                ]
            print(*fields, sep='\t')


# Get the command line arguments
if __name__ == "__main__":
    theArgs = getArgs(sys.argv[1:])

    if 'annotationMode' in theArgs and theArgs['annotationMode'] == 'vep':
        annotationMode = 'vep'
    elif 'annotationMode' in theArgs and theArgs['annotationMode'] == 'phenotype':
        annotationMode = 'phenotype'
    else:
        annotationMode = 'clinvar'

    if theArgs['regionParts'] is None:
        regionBins = parseRegionIntoBins(theArgs['region'], theArgs['binLength'], annotationMode)
    else:
        regionBins = parseRegionPartsIntoBins(theArgs['region'], theArgs['regionParts'], annotationMode)

    binIndex = +0
    data = summarizeVariants(regionBins, annotationMode)

    printResults(data, annotationMode)
