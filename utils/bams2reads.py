#!/bin/env python
import sys
import pysam
import re
import os

# for each coordinate, fetch the read names of two bam files to intersect.
# TODO: make sure pysam is zero based.

fsam1 = '/datasets/public/IlluminaPlatinum/ERP001960/BAM/bwa-mem/NA12878_1.bwa.sorted.recal.fixRG.bam'
fsam2 = '/datasets/public/IlluminaPlatinum/ERP001960/BAM/SJMAP/SJNORM016778_G1-NA12878.bam'
coords = sys.argv[1]
print_names = True
#coords2 = "chr1:12085114-12085114"

def coords_decompose(coo, genome_build):
    chrom, p0, p1 = re.split('[:-]', coo)
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    if genome_build == "hg19":
        chrom = "chr" + chrom
    return chrom, int(p0), int(p1)

# fetch return generator of alignedread reads
def fetch_read_names(sam, chrom, p0, p1):
    reads = sam.fetch(chrom, p0-1, p1)
    return [r.qname for r in reads]

def bam2qname(fsam, coo, genome_build):
    sam = pysam.Samfile(fsam, 'rb')
    chrom, p0, p1 = coords_decompose(coo, genome_build)
    return fetch_read_names(sam, chrom, p0, p1)

def main(print_names):
    name_set1 = set(bam2qname(fsam1, coords, 'hg19'))
    name_set2 = set(bam2qname(fsam2, coords, 'b37'))
    c_shared = name_set1 & name_set2
    c1 = name_set1 - name_set2
    c2 = name_set2 - name_set1
    print '\n# Private to %s at %s\nCQ:\t%d' % (os.path.basename(fsam1), coords, len(c1))
    if print_names: print 'NQ:\t' + ' '.join(sorted(c1))
    print '\n# Private to %s at %s\nCQ:\t%d' % (os.path.basename(fsam2), coords, len(c2))
    if print_names: print 'NQ:\t' + ' '.join(sorted(c2))
    print '\n# Shared:\nCQ:\t%d' % len(c_shared)
    if print_names: print 'NQ:\t' + ' '.join(sorted(c_shared))

if __name__ == '__main__':
    main(print_names)

