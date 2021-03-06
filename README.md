# bwa aln pipeline    
- bwa mapping + dedup    
- bwa hybrid mapping + dedup + split    
The pipeline works for single references such as GRCh37 and GRCh38, as well as hybrid genomes (e.g. mm9+dm3).    

## dependencies
+ bwa
+ python 2.7
+ SJM (simple job manager)
+ bgzip and tabix 1.2.1
+ LSF batch system (or other)


## set up module
export MODULEPATH=$MODULEPATH:/dir/to/bwamapping/modulefiles
module load bwamapping

## run module
run_bwa_aln.py -r1 ../fq/*_R1.fastq -r2 ../fq/*_R2.fastq -o b37.bam -O `pwd` --id NA12878 --pl ILN --sm NA12878 --lb GIAB --tmp /scratch_space -j b37.sjm

## understand ID and SM in RG tag    
ID:     
- ID is taged for each read in teh BAM     
- ID do not need to carry any meaning as it is used as uniq key.    
- Can be used to distinguish the reads from different experiments. for example if your initial sequencing
does not have enough coverage, and a topoff was done to get more reads. ID could be "NA12878" and "NA12878-topoff" to distinguish every reads.    
    

SM:     
- Every BAM is encouraged to have one SM tag, unless in rare scenario you need to merge numtiple samples together.     
- SM is used by GATK as sample name in variant call step (write to VCF file). if there are multiple SM tag in a BAM GATK HC caller (gvcf mode) would be confused.    
- SM can be identical to ID for one sample, one experiemnts cases.    

## SJM
https://github.com/StanfordBioinformatics/SJM    
