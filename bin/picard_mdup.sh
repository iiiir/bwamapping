#!/bin/bash -eu

>&2 echo "*** Mark/Remove duplicates ***"

if [ $# -lt 2 ]
then 
	>&2 echo "Usage: $0 <out.bam> <in.1.bam> [in.2.bam] [in.3.bam] ..."
	exit 1
fi

rmdup="false"

o=`cd \`dirname $1\`; pwd`/`basename $1`
shift

f=""
for bam in $@; do
	b=`cd \`dirname $bam\`; pwd`/`basename $bam`
	f="$f I=$b"
done

>&2 echo ">>> Marking duplicates"
cmd="java -Xms5g -Xmx5g -jar $PICARDPATH/picard.jar MarkDuplicates\
	TMP_DIR=$JAVATMP \
	$f \
	O=${o} \
	M=${o/.bam/.metrics} \
	VALIDATION_STRINGENCY=SILENT \
	ASSUME_SORTED=true \
	REMOVE_DUPLICATES=$rmdup"
eval $cmd

>&2 echo ">>> Indexing $o"
cmd="sam_index.sh $o"
eval $cmd

>&2 echo "*** Finished removing duplicates ***"
