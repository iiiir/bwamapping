#!/bin/bash
[[ $# -lt 2 ]] && echo "$0 <prefix> <vcf1> <vcf2>" && exit 1

vcf1=$1
[[ ! $vcf1 = *.vcf.gz ]] && echo "$vcf1 must be xxx.vcf.gz " && exit 0 
vcf1_name=`basename $vcf1 .vcf.gz`
vcf2=$2
[[ ! $vcf2 = *.vcf.gz ]] && echo "$vcf2 must be xxx.vcf.gz " && exit 0
vcf2_name=`basename $vcf2 .vcf.gz`

prefix_name="temoperary_file_name_that_you_dont_care"

cmd="vcf-isec -p $prefix_name $vcf1 $vcf2"
echo "Your command is:" $cmd
eval $cmd

mv ${prefix_name}0.vcf.gz $vcf1_name.private.vcf.gz
mv ${prefix_name}1.vcf.gz $vcf2_name.private.vcf.gz
mv ${prefix_name}0_1.vcf.gz $vcf1_name.$vcf2_name.shared.vcf.gz
mv ${prefix_name}_README $vcf1_name.$vcf2_name.README
