#!/bin/bash -eu

# 07/11/2016
# sw
# adding ATTRIBUTES_TO_RETAIN

>&2 echo "*** Merge unaligned and aligned bams ***"

if [ $# -lt 3 ]
then
    >&2 echo "Usage: $0 <out.bam> <aln.bam> [u.bam]"
    exit 1
fi

o=`cd \`dirname $1\`; pwd`/`basename $1`
shift

aln=`cd \`dirname $1\`; pwd`/`basename $1`
shift

u=`cd \`dirname $1\`; pwd`/`basename $1`
shift

# by default, tags starts with XYZ must be specified
# all bwa tags were kept
cmd="java -Xms5g -Xmx5g -jar /hpcf/apps/picard/install/2.0.1/picard.jar MergeBamAlignment \
		R=$ref_genome \
		O=$o \
		CREATE_INDEX=true \
		PAIRED_RUN=true \
		ADD_MATE_CIGAR=true \
		CLIP_ADAPTERS=false \
		CLIP_OVERLAPPING_READS=true \
		INCLUDE_SECONDARY_ALIGNMENTS=true \
		MAX_INSERTIONS_OR_DELETIONS=-1 \
		PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
		VALIDATION_STRINGENCY=LENIENT \
		ATTRIBUTES_TO_RETAIN=X0 \
		ATTRIBUTES_TO_RETAIN=X1 \
		ATTRIBUTES_TO_RETAIN=XN \
		ATTRIBUTES_TO_RETAIN=XM \
		ATTRIBUTES_TO_RETAIN=XO \
		ATTRIBUTES_TO_RETAIN=XG \
		ATTRIBUTES_TO_RETAIN=XT \
		ATTRIBUTES_TO_RETAIN=XA \
		ATTRIBUTES_TO_RETAIN=XS \
		ATTRIBUTES_TO_RETAIN=XF \
		ATTRIBUTES_TO_RETAIN=XE \
		TMP_DIR=$JAVATMP \
		UNMAPPED_BAM=$u \
		ALIGNED_BAM=$aln"
echo $cmd
eval $cmd
