"""
process bowtie2 output to make chimera.

Note. Bowtie2 does not produce supplementary alignments. It has one primary and
multiple secondary alignments ("FLAG 256"). In the unsorted sam output, the
first is always the primary alignment. If any of the secondary alignment can
combine with the primary, with at most overlap (eg. 2nt) in the query,
then we consider the pair as a chimera, and modify the two alignments with tags
'ch:A:1\tSA:Z:A'. If multiple loci for secondary alignments can be matched to
the primary, keep one for simplicity (currently STAR only reports 1 chimera for
each read). Only reads without linkers are considered now in order to be
comparable to a typical STAR run. In reality, the linkers can be easily removed
before the STAR mapping step.

For the analysis of RNA virus co-infection data, we need to examine the details
of the mapping scores and mutations. 

example command:
cd /home/zhipeng/HEK1/ #on changrila2
python ~/bin/bowtie2chim.py 

cd ~/Desktop/ #on my own Mac
python ~/Documents/scripts/duplex/bowtie2chim.py 3 \
AMT_Stress_trim_nodup_bc07_bt2aligater100000.sam \
AMT_Stress_trim_nodup_bc07_bt2aligater100000out.sam

python ~/Documents/scripts/duplex/bowtie2chim.py 3 \
AMT_Stress_trim_nodup_bc07_hg38priold1Aligned_ACTBsort.sam \
AMT_Stress_trim_nodup_bc07_hg38priold1Aligned_ACTBout.sam

"""


import sys
import itertools
import re

if len(sys.argv) < 5:
    print("Usage: python bowtie2chim.py overlap distance inputsam outputsam")
    print("overlap: overlap length of the 2 segments in the read")
    print("distance: distance of the 2 segments in the read")
    print("Not sure if the first secondary alignment is the second best")
    print("how to determine which pair of chimeric to keep?")
    sys.exit()
overlap = int(sys.argv[1])
distance = int(sys.argv[2])
samfile = open(sys.argv[3], 'r')
outfile = open(sys.argv[4], 'w')

def mergeCIGAR(CIGAR): 
    #calculate the bounds of the query, excluding the terminal SH
    #merge all operations that consume the query, i.e. MI=X
    #example: 1S2M3N4M5I6M7S -> 1S2M4M5I6M7S -> 1S17M7S (17nt query) -> [1,18]
    ops = re.findall('\d+[MISH=X]', CIGAR) #all that consume query
    newops = [ops[0]]
    for op in ops[1:]: #concatenate all internal ops that consume query [MIS=X]
        if newops[-1][-1] in "SH" and op[-1] in "MI=X": newops.append(op)
        elif newops[-1][-1]=="M" and op[-1] in "MI=X":
            newops[-1]=str(int(newops[-1][:-1])+int(op[:-1]))+"M"
        else: newops.append(op) #newops[-1][-1] == "M" and op[-1] in "SH"
    #the new ops should always contain one M with or without SH on either side
    #example: ["xM"], ["xSH","xM"], ["xM", "xSH"], ["xSH","xM", "xSH"]
    interval = []
    if len(newops)==1 or newops[0][-1]=="M": interval=[0,int(newops[0][:-1])]
    else:interval=[int(newops[0][:-1]),int(newops[0][:-1])+int(newops[1][:-1])]
    return interval
    #print("1S2M3N4M5I6M7S ->", mergeCIGAR("1S2M3N4M5I6M7S")) #output: [1, 18]


def getOverlap(a, b):
    #determine the overlap of two intervals, both exclusive
    #returns 0 if no overlap, the overlap size if overlap
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))
    #example:print(getOverlap([0,10],[10,11]))




header = ''
alignlist = [] #store all alignments for the same read
readname = ''
i = 1


for line in samfile:
    if line[0] == "@": outfile.write(line)
    else: alignlist.append(line); break

for line in samfile:
    i+=1
    if not i%1000000: print("Processed", i, "lines ...")
    if line.split()[5] == "*": continue
    readname = line.split()[0]
    if readname == alignlist[0].split()[0]: alignlist.append(line)
    else:
        #if len(alignlist) == 1: continue #ignore unique mapped
        cigar0 = alignlist[0].split()[5] #in unsorted sam, the first is primary
        intvl0 = mergeCIGAR(cigar0)
        line0 = alignlist[0].strip('\n') + '\tch:A:1\tSA:Z:A\n'
        outlist = []
        for line1 in alignlist[1:]:
            cigar1 = line1.split()[5]; intvl1 = mergeCIGAR(cigar1)
            if -distance <=max(intvl0)-min(intvl1)<=overlap \
               or -distance <=max(intvl1)-min(intvl0)<=overlap \
               and len(re.findall("M",cigar0))==1 \
               and len(re.findall("M",cigar1))==1:
                outlist.append(line0+line1.strip('\n')+'\tch:A:1\tSA:Z:A\n')
        if outlist: outfile.write(outlist[0]) #take only one pair for each read
        alignlist = [line]
samfile.close()
outfile.close()










