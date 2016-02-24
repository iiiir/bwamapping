#!/bin/env python
import sys
import numpy as np
import pysam
import matplotlib.pyplot as plt


fnsam = sys.argv[1]
fsam=pysam.Samfile(fnsam)
#N = 75891268
#NM_scores = np.arange(N)
def find_nm(lis):
    nm_val = '255'
    for k,v in lis:
        if k == "NM": nm_val = v
    return nm_val

for i,read in enumerate(fsam):
#    NM_scores[i] = find_nm(read.tags)
    print find_nm(read.tags)
#uniq_nm, inds = np.unique(NM_scores, return_index=True)

#plt.hist(NM_scores, bin=100)
#plt.show()

