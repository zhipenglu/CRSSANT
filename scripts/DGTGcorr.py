"""
DGTGcorr.py

This script calculates Pearson correlation between TG coverage and DG coverage
Here DG coverage is defined as the square-root of the product of the alignment
nummbers in the 2 DGs that support the TG. 

input file: 
Column1: TG info. Contains the original DG ID, confidence and alignment numbers
Column2: TG coverage. 
E.g.: TG:Z:28S,28S,3,0.071_1621_28S,28S,20,0.283_6441	3865

Example test on the 28S rRNA DGs and TGs (top 242 TGs):
cd ~/Desktop/
python ~/Documents/scripts/crssant/DGTGcorr.py parisTG.txt
python ~/Documents/scripts/crssant/DGTGcorr.py parisShuffle.txt
"""

import sys
import math
import scipy.stats as stats
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) <2:
    print("Usage: python DGTGcorr.py inputfile")
    sys.exit()
infile = open(sys.argv[1], 'r')
    
##### extract numbers and calculate correlation
DGcov = [] #geometric means of read numbers in the two DGs that make up each TG
TGcov = [] #read numbers for each TG.
for line in infile:
    TG1, TG2 = tuple(line.split())
    DGcov.append(math.log10((int(TG1.split('_')[1])*int(TG1.split('_')[3]))**0.5))
    TGcov.append(math.log10(int(TG2)))

##### using all TGs or the top 242?
DGcov = DGcov[:242]
TGcov = TGcov[:242]
print("Pearson correlation:", stats.pearsonr(DGcov, TGcov))

##### plot the numbers in a scatter
fig, ax = plt.subplots(figsize=(10,10))
plt.scatter(DGcov, TGcov)
m, b = np.polyfit(DGcov, TGcov, 1)
print("Correlation parameters:", m,b)
plt.plot(np.array(DGcov), m*np.array(DGcov)+b)
plt.xlim(0,4)
plt.ylim(0,4)
#further plot the correlation 
plt.savefig(sys.argv[1]+".pdf")












