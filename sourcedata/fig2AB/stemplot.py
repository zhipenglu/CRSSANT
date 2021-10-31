"""
stemplot.py. Copyright: Zhipeng Lu, 2020, zhipengluchina@gmail.com
This script calculates the size distribution of the RNA stems (each side counted
separately).


example plot:
cat hs_lsu_rrna_helix.bed hs_ssu_rrna_helix.bed > hs_rrna_helix.bed

cd /Users/lu/Documents/lulab/rrna/hsrRNA_helix
python ~/Documents/scripts/pyplots/stemplot.py hs_rrna_helix.bed \
hs_rrna_helix_dist.pdf

cd /Users/lu/Documents/chang/psoralen/examples/snRNA/hssnRNA_bphelix
python ~/Documents/scripts/pyplots/stemplot.py hs_snRNA_all_helix.bed \
hs_snRNA_all_helix_dist.pdf
"""

import matplotlib.pyplot as plt
import sys, re, matplotlib, matplotlib.pyplot as plt

if len(sys.argv)< 3:
    print("Usage: python stemplot.py inputfile outputfile")
    print("input is bed12 format, output is the pdf figure")
    sys.exit()

inputfile = open(sys.argv[1], 'r')
sizelist = [] #create the size list from the input bed file. 
for line in inputfile: sizelist+=[int(i) for i in line.split()[10].split(',')]
print("Number of stem arms:", len(sizelist))
print("Fraction of arms <=20nt:",
      len([i for i in sizelist if i <=20])/len(sizelist))
print("Fraction of arms <=10nt:",
      len([i for i in sizelist if i <=10])/len(sizelist))
print("Fraction of arms <=5nt:",
      len([i for i in sizelist if i <=5])/len(sizelist))

fig, ax = plt.subplots(figsize=(2.5,1.5))
plt.subplots_adjust(bottom=0.32)
plt.subplots_adjust(top=0.95)
plt.subplots_adjust(left=0.3)
plt.subplots_adjust(right=0.95)
n,bins,patches=plt.hist(sizelist, [x+0.5 for x in range(0,100)],color='b',
                        histtype='step',cumulative=False,normed=1,fill='b') 
#plt.setp(patches, 'facecolor', 'b', 'alpha', 0.75)
plt.xlim(-3, 40)
#plt.ylim(0, 1.1)
#plt.locator_params(nbins=5, axis='y')

max_yticks = 5
yloc = plt.MaxNLocator(max_yticks)
ax.yaxis.set_major_locator(yloc)

ax.set_xlabel("stem arm length (nt)") #, fontsize=15
ax.set_ylabel("frequency") #, fontsize=15
plt.savefig(sys.argv[2])
plt.show()
           


