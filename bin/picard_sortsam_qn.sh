#!/bin/bash -eu

>&2 echo "*** sort bam by query names ***"

if [ $# -lt 1 ]
then
    >&2 echo "Usage: $0 <aln.bam> [output.bam]"
    exit 1
fi

bam=`cd \`dirname $1\`; pwd`/`basename $1`; shift
[[ $# -gt 0 ]] && o=$1 || o=`cd \`dirname $bam\`; pwd`/`basename $bam .bam`.qn.bam

[[ -f $o ]] && echo "$o existed" && exit 0

cmd="java -jar /hpcf/apps/picard/install/2.0.1/picard.jar SortSam \
		TMP_DIR=$JAVATMP \
		I=$bam \
		O=$o \
		SORT_ORDER=queryname"
echo $cmd
eval $cmd
