#!/bin/bash -eu

echo "*** Sorting BAM by position ***"

if [ $# -lt 1 ]
then 
	echo "Usage: $0 <bam> [memory in GB]"
	exit 1
fi

f=`cd \`dirname $1\`; pwd`/`basename $1`
o=${f/.bam/.sorted}

gmem=5G
if [ $# -gt 2 ]
then
	gmem=$3
fi

echo ">>> Sorting on BAM $f"
samtools sort -m $gmem $f $o && samtools index $o.bam

echo "*** Finished sorting BAM by position ***"
