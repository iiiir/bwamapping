#!/bin/bash

[[ $# -lt 2 ]] && echo "$0 <to_be_liftover.vcf> <to_hg19|to_b37>" && exit 0

ivcf=`realpath.sh $1`
oform=$2
ovcf=`basename ${ivcf/.vcf/.${oform#to_}.vcf}`

if [ $oform = "to_hg19" ]; then
	#From:
	ref_genome="/datasets/public/genomes/hsapiens/hg19/SNPS/gatk_bundle/ref_genome/GRCh37.fa"
	CHAINFILE="/datasets/public/genomes/hsapiens/hg19/SNPS/gatk_bundle/liftover/b37tohg19.chain"
	# to
	to_ref="/datasets/public/genomes/hsapiens/hg19/SNPS/gatk_bundle/ref_genome/hg19.fa"
	DICTFILE="/datasets/public/genomes/hsapiens/hg19/SNPS/gatk_bundle/ref_genome/hg19.dict"
elif [ $oform = "to_b37" ]; then
	#From:
    ref_genome="/datasets/public/genomes/hsapiens/hg19/SNPS/gatk_bundle/ref_genome/hg19.fa"
    CHAINFILE="/datasets/public/genomes/hsapiens/hg19/SNPS/gatk_bundle/liftover/hg19tob37.chain"
    # to
	to_ref="/datasets/public/genomes/hsapiens/hg19/SNPS/gatk_bundle/ref_genome/GRCh37.fa"
    DICTFILE="/datasets/public/genomes/hsapiens/hg19/SNPS/gatk_bundle/ref_genome/GRCh37.dict"
else
	echo "$0 <to_be_liftover.vcf|to_be_liftover.vcf.gz> <to_hg19|to_b37>"
	exit 0
fi

>&2 echo "*** LiftOver VCF file $oform ***"
>&2 echo ">>> Step 1: LiftOver variants $oform"
cmd="java -Xms3g -Xmx3g -XX:ParallelGCThreads=4 -Djava.io.tmpdir=$JAVATMP \
	-jar ${GATKPATH}/GenomeAnalysisTK.jar \
	-T LiftoverVariants \
	-R $ref_genome \
	-chain $CHAINFILE \
	-dict $DICTFILE \
	-V $ivcf \
	-o $ovcf.temp"

echo $cmd

>&2 echo ">>> Step 2: Filter variants to correct header"
cmd="java -Xms3g -Xmx3g -XX:ParallelGCThreads=4 -Djava.io.tmpdir=$JAVATMP \
	-jar ${GATKPATH}/GenomeAnalysisTK.jar \
	-T FilterLiftedVariants \
	-R $to_ref \
	-V $ovcf.temp \
	-o $ovcf"
echo $cmd

>&2 echo "*** Finished liftOver ***"

