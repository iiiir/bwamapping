#!/bin/bash -e

>&2 echo "*** sort bam by query names ***"

if [ $# -lt 1 ]
then
    >&2 echo "Usage: $0 <in.u.bam> [out.bam]"
    exit 1
fi

inbam=$1; shift
[[ $# -gt 0 ]] && obam=$1 || obam=`basename $inbam .u.bam`.bam

cmd="java -jar /hpcf/apps/picard/install/2.0.1/picard.jar SortSam \
		TMP_DIR=$JAVATMP \
		I=$inbam \
		O=$obam \
		SORT_ORDER=queryname"
echo $cmd
eval $cmd

# qc
read_count1=`samtools view $inbam | wc -l`
read_count2=`samtools view $obam | wc -l`
if [[ ! $read_count1 = $read_count2 ]]; then
	>&2 echo "$inbam -> $obam"
	>&2 echo "$read_count1 reads -> $read_count2 reads"
	exit -1
fi
>&2 echo "*** finished sort bam by query names ***"
