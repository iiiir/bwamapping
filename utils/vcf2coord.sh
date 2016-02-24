#!/bin/bash
[[ ! $# -eq 1 ]] && echo "$0 1.vcf" && exit 1
f=$1
awk -F"\t" '($1 !~ /^#/) {print $1":"$2"-"$2}' $f
