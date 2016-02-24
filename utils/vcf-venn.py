#!/bin/env python
from matplotlib import pyplot as plt
from matplotlib_venn import venn3, venn3_circles
import sys
import subprocess
import os

vcftools_path = "/hpcf/apps/vcftools/install/0.1.12b/bin"

def parse_vcf_compare(fn):
    '''Open summary file from vcf-compare and parse numbers'''
    #lis = [l.split('\t')[1] for l in open(fn) if l.startswith('VN')]
    n101, n001, n011, n010, n100, n110, n111 = [int(l.split('\t')[1]) for l in open(fn) if l.startswith('VN')]
    return n101, n001, n011, n010, n100, n110, n111

def parse_vcf_compare_percentage(fn):
    n101, n001, n011, n010, n100, n110, n111 = [l.strip().split(' ')[-1][1:-1] for l in open(fn) if l.startswith('VN')]
    return n001, n101, n011, n111


def vcf_compare(vcfs):
    cmd = [os.path.join(vcftools_path,'vcf-compare')] + [os.path.realpath(v) for v in vcfs] + ['>'] + ['mytest_vcf_compare.summary']
    subprocess.Popen(cmd) 

def main(fn):
    plt.figure(figsize=(8.5,11))
    n101, n001, n011, n010, n100, n110, n111 = parse_vcf_compare(fn)
    n001s, n101s, n011s, n111s = parse_vcf_compare_percentage(fn)
    l001, l101, l011, l111 = ['%d(%s)'%(a,b) for a,b in zip((n001, n101, n011, n111), (n001s, n101s, n011s, n111s))]
    mysets = (n100,n010,n110,n001,n101,n011,n111)
    mylables = ('FreeBayes', 'GATiK (UG)', 'Plantinum v0.7')
    v = venn3(subsets    = mysets, 
              set_labels = mylables)
    #venn3_circles(subsets=mysets, linestyle='dashed') # solid
    plt.title("Venn diagram")
    #v.get_patch_by_id('100').set_color('white')
    #v.get_patch_by_id('010').set_color('white')
    #v.get_patch_by_id('110').set_color('white')
    v.get_label_by_id('001').set_text(l001)
    v.get_label_by_id('011').set_text(l011)
    v.get_label_by_id('101').set_text(l101)
    v.get_label_by_id('111').set_text(l111)
    plt.show()

if __name__ == "__main__":
    fn = sys.argv[1]
    main(fn)
