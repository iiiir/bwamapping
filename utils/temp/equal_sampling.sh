#!/bin/bash

myfile=$1

while read line; do
    bams2reads.py $line | grep CQ | cut -f2 | paste -sd " " >> $myfile.out
done<$myfile
    
