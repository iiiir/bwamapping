#!/bin/bash -eu

echo "*** Aligning reads using BWA aln/PE algorithm ***"
echo ">>> Caveat: pipe might broken but file will be written!"
echo ">>> e.g. get empty out.bam because fastq file does not exist"

if [ $# -lt 1 ]
then 
	echo "Usage: $0 <out.bam> <u.bam> <u.bam> [RG tag]"
	exit 1
fi

o=`cd \`dirname $1\`; pwd`/`basename $1`; shift
[[ $o = *.bam ]] && o=${o%.bam}

q1=`cd \`dirname $1\`; pwd`/`basename $1`
q2=`cd \`dirname $2\`; pwd`/`basename $2`

[[ ! -f $q1 ]] && echo ">>> Does not exist: $q1" && exit 1
[[ ! -f $q2 ]] && echo ">>> Does not exist: $q2" && exit 1

## update 07/23/2015: change -R to -r for RG tag
optRG=""
lastArg=${BASH_ARGV[0]}
if [[ $lastArg =~ "@RG" ]]
then
        optRG="-r $lastArg"
fi

cmd="bwa sampe $ref_genome \
            <(bwa aln -1 -t 4 $ref_genome $q1) \
            <(bwa aln -2 -t 4 $ref_genome $q1) $q1 $q1 $optRG | \
            samtools view -b - | \
            samtools sort - $o"
echo $cmd
eval $cmd

[[ -f $o.bam ]] && samtools index $o.bam || (echo ">>> not found: $o.bam"; exit 1)

echo "*** Finished aligning reads ***"
