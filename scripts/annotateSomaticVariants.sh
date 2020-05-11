#!/bin/bash
#SJG Jan2020

vcfUrl=$1
selectedSamples=$2
regions=$3
somaticFilterPhrase=$4
genomeBuildName=$5
vepCacheDir=$6

runDir=$PWD
tempDir=$(mktemp -d)
cd $tempDir

echo -e "$regions" >> regions.txt

<<'OLD'
echo -e "$normalSampleIdxs" >> normals.txt
echo -e "$tumorSampleIdxs" >> tumors.txt

qualPhrase="QUAL>${qualCutoff}"
depthPhrase="INFO/DP>${depthCutoff}"

cmpdNormalPhrase="("
n_count=0
while read idx; do
	currPhrase="("
	normalCountPhrase="FORMAT/AO[${idx}:0]<=${normalCountCutoff}"
	normalFreqPhrase="(FORMAT/AO[${idx}:0]/(FORMAT/AO[${idx}:0]+FORMAT/RO[${idx}:0]))<=${normalAfCutoff}"
	if (("$n_count" > 0)); then
		currPhrase="|("
	fi
	currPhrase="${currPhrase}${normalCountPhrase}|${normalFreqPhrase})"
	cmpdNormalPhrase="${cmpdNormalPhrase}${currPhrase}"
	((n_count=n_count+1));
done < normals.txt
cmpdNormalPhrase="${cmpdNormalPhrase})"

cmpdTumorPhrase="("
t_count=0
while read idx; do
        currPhrase="("
        tumorCountPhrase="FORMAT/AO[${idx}:0]>=${tumorCountCutoff}"
        tumorFreqPhrase="FORMAT/AO[${idx}:0]/(FORMAT/AO[${idx}:0]+FORMAT/RO[${idx}:0])>=${tumorAfCutoff}"
	if (("$t_count" > 0)); then
                currPhrase="|("
        fi
        currPhrase="${currPhrase}${tumorCountPhrase}|${tumorFreqPhrase})"
        cmpdTumorPhrase="${cmpdtumorPhrase}${currPhrase}"
	((t_count=t_count+1))
done < tumors.txt
cmpdTumorPhrase="${cmpdTumorPhrase})"


normalCountPhrase="AC[${normalSampleIdx}]<=${normalCountCutoff}"
normalAfPhrase="(AF[${normalSampleIdx}]<=${normalAfCutoff}||(AC[${normalSampleIdx}]/AN)<=${normalAfCutoff})"

#format tumor query pieces
tumorCountPhrase="("
tumorAfPhrase="("
for ((i=0; i<totalSampleNum; i++)); do
	if ((i != normalSampleIdx)); then
		tumorCountPhrase="${tumorCountPhrase}AC[${i}]>=${tumorCountCutoff}"
		tumorAfPhrase="${tumorAfPhrase}((AF[${i}]>=${tumorAfCutoff})||((AC[${i}]/AN)<=${tumorAfCutoff}))"
		if ((i < totalSampleNum-1)); then
4		tumorCountPhrase="${tumorCountPhrase}||"
			tumorAfPhrase="$tumorAfPhrase||"
		fi
	fi
done
tumorCountPhrase="${tumorCountPhrase})"
tumorAfPhrase="${tumorAfPhrase})"

OLD

#format final query
#queryPhrase="${qualFilterPhrase}&${normalFilterPhrase}&${tumorFilterPhrase}"

somaticFilterPhrase="FORMAT/AO[0:0]<2&FORMAT/AO[1:0]>5|FORMAT/AO[2:0]>5|FORMAT/AO[3:0]>5"

echo "${somaticFilterPhrase}"

fastaPath="/home/ubuntu/data/references/GRCh37/human_g1k_v37_decoy_phix.fasta"
if [ "${genomeBuildName}" = "GRCh38" ]; then
	fastaPath="/home/ubuntu/data/references/GRCh38/human_g1k_v38_decoy_phix.fasta"
fi

vepArgs="--assembly $genomeBuildName --format vcf --allele_number --dir_cache $vepCacheDir"

bcftools view -s $selectedSamples $vcfUrl | \
    bcftools filter -t $regions - | \
    bcftools norm -m - -w 10000 -f $fastaPath - | \
    bcftools filter -i $somaticFilterPhrase
    #vep $vepArgs

rm -rf $tempDir
cd $runDir

