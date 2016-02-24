#!/bin/bash -eu

echo "*** Aligning reads using BWA aln/PE algorithm ***"
echo ">>> Caveat: pipe might broken but file will be written!"
echo ">>> e.g. get empty out.bam because fastq file does not exist"

if [ $# -lt 1 ]
then 
	echo "Usage: $0 <out.bam> <q1> <q2> <ref.fa> [RG tag]"
	exit 1
fi

# define output
o=`cd \`dirname $1\`; pwd`/`basename $1`; shift
[[ $o = *.bam ]] && o=${o%.bam}

# define input reads
q1=`cd \`dirname $1\`; pwd`/`basename $1`; shift
q2=`cd \`dirname $1\`; pwd`/`basename $1`; shift
if [[ $q1 = *.bam ]]; then
	input_format="ubam"
	sequencing=$q2 # if ubam provided, must specify SE or PE
fi

[[ ! -f $q1 ]] && echo ">>> Does not exist: $q1" && exit 1

if [[ ( $q1 = *.fastq ) || ( $q1 = *.fq ) || ( $q1 = *.gz ) ]]; then
	input_format="fastq"
	if [[ $q1 = $q2 ]]; then
		sequencing="SE"
	else
		sequencing="PE"
		[[ ! -f $q2 ]] && echo ">>> Does not exist: $q2" && exit 1
	fi
fi

# define reference
ref_fasta=`cd \`dirname $1\`; pwd`/`basename $1`; shift

optRG=""
lastArg=${BASH_ARGV[0]}
if [[ $lastArg =~ "@RG" ]]
then
        optRG="-r $lastArg"
fi

if [[ $input_format = "ubam" ]]; then
	echo ">>> Your input is ubam"
	if [[ $sequencing = "SE" ]]; then
		echo ">>> Your sequencing type is single-end"
		cmd="bwa aln -b -0 -t 4 $ref_fasta $q1 | \
			 bwa samse $ref_fasta - $q1 $optRG | \
			 samtools view -b - | \
			 samtools sort - $o"
	else
		echo ">>> Your sequencing type is pair-end"
		cmd="bwa sampe $ref_fasta \
			<(bwa aln -b -1 -t 4 $ref_fasta $q1) \
			<(bwa aln -b -2 -t 4 $ref_fasta $q1) $q1 $q1 | \
			samtools view -b - | \
			samtools sort - $o"		
	fi
else
	echo ">>> Your input is fastq"
	if [[ $sequencing = "SE" ]]; then
		echo ">>> Your sequencing type is single-end"
		cmd="bwa aln -t 4 $ref_fasta $q1 | \
			bwa samse $ref_fasta - $q1 $optRG | \
			samtools view -b - | \
			samtools sort - $o"
	else
		echo ">>> Your sequencing type is pair-end"
		cmd="bwa sampe $ref_fasta \
			<(bwa aln -t 4 $ref_fasta $q1) \
			<(bwa aln -t 4 $ref_fasta $q2) $q1 $q2 | \
			samtools view -b - | \
			samtools sort - $o"
	fi
fi
echo $cmd
eval $cmd

[[ -f $o.bam ]] && samtools index $o.bam || (echo ">>> not found: $o.bam"; exit 1)

echo "*** Finished aligning reads ***"
