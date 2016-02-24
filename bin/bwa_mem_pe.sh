#!/bin/bash -eu

echo "*** Aligning reads using BWA MEM algorithm ***"

if [ $# -lt 4 ]
then 
	echo "Usage: $0 <out_bam_prefix> <fastq1> <fastq2> <RG tag>"
	exit 1
fi

o=`cd \`dirname $1\`; pwd`/`basename $1`; shift

q1=`cd \`dirname $1\`; pwd`/`basename $1`
q2=`cd \`dirname $2\`; pwd`/`basename $2`

[[ ! -f $q1 ]] && echo ">>> Does not exist: $q1" && exit 1
[[ ! -f $q2 ]] && echo ">>> Does not exist: $q2" && exit 1

optRG=""
lastArg=${BASH_ARGV[0]}
if [[ $lastArg =~ "@RG" ]]
then
        optRG="-R $lastArg"
fi

cmd="bwa mem -M $ref_genome $q1 $q2 -t 4 $optRG | samtools view -b - | samtools sort - $o"
eval $cmd

sam_index.sh $o.bam

echo "*** Finished aligning reads ***"
