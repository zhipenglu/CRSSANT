"""
softreverse.py
Zhipeng Lu, zhipengluchina@gmail.com

STAR mapping is biased against backward gapped compared to normal chimeras. 
This script reverses the softclipped alignments in cont.sam for another round
of STAR mapping to rescue reads that are poorly mapped in the first STAR round.

In theory, there should be <= 2 softclips. Here we move both as long as the
longer one is longer than minlen (e.g. 5)

example:
/project/zhipengl_72/zhipengl/HEK1
AMT_Stress_trim_nodup_bc07_hg38priAligned_sorted_ACTB_include.sam

Input sam file containing all alignments.
convert only the ones with a softclip and no secondary alignments. These are
likely to contain backward alignments that are

e.g. extract alignments mapped to ACTB and check them directly. like Fig. 2GH.

example command:
cd /Users/lu/Documents/lulab/projects/computational/starmapping
python ~/Documents/scripts/crssant/github/softreverse.py \
AMT_Stress_trim_nodup_bc07_hg38priAligned_sorted_ACTB_include.sam \
AMT_Stress_trim_nodup_bc07_hg38priAligned_sorted_ACTB_include_softrev.sam

only take primary alignments. 

copy the function from a different file
For the simulation I need to make both forward and backward alignments. We do
not need to make a specific

What is the best approach? The ACTB region


"""

import sys
import re



if len(sys.argv) < 3:
    print("Usage: python softreverse.py insam outfastq")
    sys.exit()

minlen = 5 #minimal length of the softclipped shorter segment
insam = open(sys.argv[1], 'r')
outfastq = open(sys.argv[2], 'w')

for line in insam:
    if line[0] == '@': continue
    record = line.split()
    FLAG = int(record[1])
    if FLAG>=256 and bin(FLAG)[-9]=='1' or \
       len(record)>=21 and "SA:Z" in line.split()[20]:
        continue #ignore secondary (FLAG=256) and chiastic alignments (SA:Z). 
    CIGAR = line.split()[5]
    subMS = re.findall('\d+[MS]', CIGAR) #substrings for M and S
    softs = re.findall('\d+S', CIGAR)
    softslen = [int(i[:-1]) for i in softs]
    if softslen and max(softslen) >= minlen:
        record = line.split()
        QNAME, SEQ, QUAL = record[0], record[9], record[10]
        #print(CIGAR, SEQ)
        if 'S' in subMS[0] and 'S' not in subMS[-1]:
            SEQ = SEQ[softslen[0]:] + SEQ[:softslen[0]]
            QUAL = QUAL[softslen[0]:] + QUAL[:softslen[0]]
        elif 'S' not in subMS[0] and 'S' in subMS[-1]:
            SEQ = SEQ[-softslen[0]:] + SEQ[:-softslen[0]]
            QUAL = SEQ[-softslen[0]:] + SEQ[:-softslen[0]]
        elif 'S' in subMS[0] and 'S' in subMS[-1]:
            SEQ = SEQ[-softslen[1]:] + \
                  SEQ[softslen[0]:-softslen[1]] + SEQ[:softslen[0]] 
            QUAL = QUAL[-softslen[1]:] + \
                   QUAL[softslen[0]:-softslen[1]] + QUAL[:softslen[0]]             
        #print(CIGAR, SEQ)
        outfastq.write('@'+QNAME + '\n' + SEQ + '\n+\n' + QUAL + '\n')

insam.close()
outfastq.close()










