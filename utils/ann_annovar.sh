#!/bin/bash

inbed=`realpath $1`
outann=${inbed/.BED/.ann}
convert2annovar.pl -format bed $inbed > $outann
annotate_variation.pl --geneanno -buildver mm9 $outann mousedb/
