#!/bin/bash -eu

echo "*** Aligning reads using BWA aln/PE algorithm ***"
echo ">>> Caveat: pipe might broken but file will be written!"
echo ">>> e.g. get empty out.bam because fastq file does not exist"

if [ $# -lt 1 ]
then 
	echo "Usage: $0 <out.bam> <u.bam> <u.bam>"
	exit 1
fi

obam=`cd \`dirname $1\`; pwd`/`basename $1`; shift

q1=`cd \`dirname $1\`; pwd`/`basename $1`
q2=`cd \`dirname $2\`; pwd`/`basename $2`
[[ ! -f $q1 ]] && echo ">>> Does not exist: $q1" && exit 1
[[ ! -f $q2 ]] && echo ">>> Does not exist: $q2" && exit 1

cmd="bwa sampe $ref_genome \
            <(bwa aln -b -1 -t 4 $ref_genome $q1) \
            <(bwa aln -b -2 -t 4 $ref_genome $q1) $q1 $q1 | \
            samtools view -b - | \
            java -Xmx10g -jar $PICARDPATH/picard.jar SortSam \
				TMP_DIR=$JAVATMP \
				SORT_ORDER=queryname \
				VALIDATION_STRINGENCY=LENIENT \
				I=/dev/stdin \
				O=$obam"
echo $cmd
eval $cmd

echo "*** Finished aligning/query_name sorting reads ***"
