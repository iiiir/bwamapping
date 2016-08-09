#!/bin/bash -eu

>&2 echo "*** Merge unaligned and aligned bams ***"

if [ $# -lt 1 ]
then
    >&2 echo "Usage: $0 <aln.bam> [query_name]"
    exit 1
fi

bam=`cd \`dirname $1\`; pwd`/`basename $1`
shift
o=`basename $bam .bam`

[[ $# -gt 0 ]] && sort_order="queryname" || sort_order="coordinate"

cmd="java -Xms8g -Xmx8g -jar /hpcf/apps/picard/install/2.0.1/picard.jar SortSam \
		TMP_DIR=$JAVATMP \
		I=$bam \
		O=$o.sort.bam \
		SORT_ORDER=$sort_order"
echo $cmd
eval $cmd
