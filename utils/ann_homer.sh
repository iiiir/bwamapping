#!/bin/bash

bedin=$1
bedstats=`basename ${bedin%.*}.stats`
bedout=`basename ${bedin%.*}.tsv`
genome=hg19
[[ ! -z $2 ]] && genome=$2

# run command
annotatePeaks.pl $bedin $genome -annStats $bedstats > $bedout 2> $bedstats.txt
