#!/bin/bash -eu

echo "*** Aligning reads using BWA MEM algorithm ***"

if [ $# -lt 2 ]
then 
	echo "Usage: $0 <o.bam> <u.bam>"
	exit 1
fi

o=`cd \`dirname $1\`; pwd`/`basename $1`; shift
ubam=`cd \`dirname $1\`; pwd`/`basename $1`

[[ ! -f $ubam ]] && echo ">>> Does not exist: $ubam" && exit 1

cmd1="java -Xms5g -Xmx5g -jar $PICARDPATH/picard.jar SamToFastq \
	INPUT=$ubam \
	FASTQ=/dev/stdout \
	INTERLEAVE=true | \
	bwa mem -p -M $ref_genome /dev/stdin | \
	samtools view -b - > $o.mem.bam"
eval $cmd1

# default: sort by coordinate and index
cmd2="java -Xms5g -Xmx5g -jar $PICARDPATH/picard.jar MergeBamAlignment \
        R=$ref_genome \
        UNMAPPED_BAM=$ubam \
        ALIGNED_BAM=$o.mem.bam \
        O=$o \
        TMP_DIR=$JAVATMP \
        MAX_INSERTIONS_OR_DELETIONS=-1 \
        VALIDATION_STRINGENCY=LENIENT"
eval $cmd2

rm $o.mem.bam

echo "*** Finished aligning reads ***"
