#!/bin/env python3
import sys, os

class VCF():
    def __init__(self,tsv_row):
        ''' TSV_ROW is a list from csv reader module'''
        # required fields:
        # CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO
        self.GENE, self._QUAL, self.SAMPLE, self.CHROM, self.POS, \
        self.CLASS, self.AA, self.PGI, self.RNA,self.TUMER_AC, \
        self.TUMOR_DP,self._NORM_AC, self._NORM_DP, self.REF, self.ALT, \
        self.FLANK, self.DUP, self.STATE = tsv_row

    @property
    def QUAL(self):
        return '30' if self._QUAL == 'SJHQ' else '.'
    @QUAL.setter
    def QUAL(self, value):
        self._QUAL = value

    @property
    def NORM_AC(self):
        return int(self._NORM_AC)
    
    @property
    def NORM_DP(self):
        return int(self._NORM_DP)

    @property
    def DP(self):
        return self.NORM_DP

    @property
    def NORM_AF(self):
        return self.NORM_AC/self.NORM_DP 

    @property
    def AC(self):
        if self.NORM_AF == 1:
            return 2
        elif self.NORM_AF == 0:
            return 0
        else:
            return 1

    @property
    def AN(self):
        return 2

    @property
    def AF(self):
        return self.AC / self.AN

    @property
    def ID(self):
        return "."

    @property
    def FILTER(self):
        return 'PASS' if self._QUAL == 'SJHQ' else "."

    @property
    def INFO(self):
        return 'DP=%d' % (self.DP)
    
    @property
    def FORMAT(self):
        return "GT:AD:DP"

    @property
    def GENO(self):
        if self.AF == 1:
            return '1/1'
        elif self.AF == 0:
            return '0/0'
        else:
            return '0/1'

    @property
    def GENOTYPE(self):
        return '%s:%d,%d:%d' % (self.GENO, self.NORM_DP - self.NORM_AC, self.NORM_AC, self.NORM_DP)    

    def __repr__(self):
        return '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (self.CHROM,self.POS,self.ID,self.REF,self.ALT,self.QUAL,self.FILTER,self.INFO, self.FORMAT, self.GENOTYPE)

    def __str__(self):
        return '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (self.CHROM,self.POS,self.ID,self.REF,self.ALT,self.QUAL,self.FILTER,self.INFO, self.FORMAT, self.GENOTYPE)

def vcf_header(sample):
    return '''\
##fileformat=VCFv4.1
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
##phasing=partial
##INFO=<ID=AC,Number=1,Type=Integer,Description="Allele count in genotypes. for each ALT allele">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FILTER=<ID=q10,Description="Quality below 10">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for Ref and ALT">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s
''' % sample

def st_variants(tsv_file, header_line=True):
    import csvs
    csvfile=open(tsv_file)
    csvreader = csv.reader(csvfile, delimiter='\t')
    if header_line:
        next(csvreader)
    return csvreader

def col_header():
    ''' Headers defined for each column'''
    pass

def sj2vcf(tsv_file):
    sample=os.path.basename(sys.argv[1])
    print(vcf_header(sample), end='')
    for row in st_variants(tsv_file):
        print(VCF(row))

def test():
    '''
    >>> sj2vcf(lin)

    '''
    #sj2vcf('/nfs_exports/genomes/1/PCGP/BucketIntermediate/SJETV/SnpDetect/test_SJETV195/germlineTest/allgenes/SJETV195_novel_germline.txt')
    sj2vcf(sys.argv[1])

def main():
    test()

if __name__ == '__main__':
    main()
