#!/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import sys

fn = sys.argv[1] #"0001_other_only_coords.reads.compare.txt" three columns of numbers
data=np.loadtxt(fn)
fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(6,6))
axes.boxplot(data)

ymin,ymax= np.percentile(data[:,2],[0,75])
print 'Median of  bwa-mem : %d' % np.percentile(data[:,0],50)
print 'Median of    bwa   : %d' % np.percentile(data[:,1],50)
print 'Median of Identical: %d' % np.percentile(data[:,2],50)
axes.set_ylim([ymin,ymax+1])
axes.set_xticklabels(['bwa-mem', 'bwa', 'Identical'])
plt.show()
