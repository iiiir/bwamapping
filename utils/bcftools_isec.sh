#!/bin/bash
[[ $# -lt 2 ]] && echo "$0 <vcf1> <vcf2>" && exit 1

vcf1=$1
[[ ! $vcf1 = *.vcf.gz ]] && echo "$vcf1 must be xxx.vcf.gz " && exit 0 
vcf1_name=`basename $vcf1 .vcf.gz`
vcf2=$2
[[ ! $vcf2 = *.vcf.gz ]] && echo "$vcf2 must be xxx.vcf.gz " && exit 0
vcf2_name=`basename $vcf2 .vcf.gz`

cmd="bcftools isec -p `pwd` $vcf1 $vcf2"
echo "Your command is:" $cmd
eval $cmd

mv 0000.vcf $vcf1_name.private.vcf
mv 0001.vcf $vcf2_name.private.vcf
mv 0002.vcf $vcf1_name.shared.vcf
mv 0003.vcf $vcf2_name.shared.vcf
