#!/usr/bin/python3
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

Jessica Severin: July 2023, modified this code to use the pysam toolkit to 
  allow reading of both BAM and SAM files to avoid the BAM->SAM conversion 
  step in the workflow.  https://pysam.readthedocs.io/en/latest/api.html

"""

import sys
import re
import pysam


if len(sys.argv) < 3:
    print("Usage: python3 softreverse.py insam outfastq")
    sys.exit()

infile = sys.argv[1]
outfile = sys.argv[2]
insam = None

if(infile.endswith("sam")):
    print(infile," is SAM")
    insam = pysam.AlignmentFile(infile, "r", require_index=False)
if(infile.endswith(".bam")):
    print(infile, " is BAM")
    insam = pysam.AlignmentFile(infile, "rb", require_index=False)

if(insam is None):
    print("Usage: python3 softreverse.py insam outfastq : must use .sam or .bam input file")
    sys.exit()

minlen = 5 #minimal length of the softclipped shorter segment
outfastq = open(sys.argv[2], 'w')

cond1_count = 0
cond2_count = 0
cond3_count = 0
total_count = 0

for align in insam:
    total_count+=1
    line = align.to_string()
    record = line.split()
    FLAG = align.flag
    #if (FLAG>=256 and bin(FLAG)[-9]=='1') or (len(record)>=21 and "SA:Z" in record[20]):
    if (FLAG>=256 and bin(FLAG)[-9]=='1') or (align.has_tag('SA')):
        continue #ignore secondary (FLAG=256) and chiastic alignments (SA:Z). 
    CIGAR = align.cigarstring
    subMS = re.findall('\d+[MS]', CIGAR) #substrings for M and S
    softs = re.findall('\d+S', CIGAR)
    softslen = [int(i[:-1]) for i in softs]
    if softslen and max(softslen) >= minlen:
        QNAME = align.query_name
        SEQ   = align.query_sequence
        QUAL  = align.qual
        #print(CIGAR, SEQ)
        if 'S' in subMS[0] and 'S' not in subMS[-1]:
            cond1_count+=1
            SEQ = SEQ[softslen[0]:] + SEQ[:softslen[0]]
            QUAL = QUAL[softslen[0]:] + QUAL[:softslen[0]]
        elif 'S' not in subMS[0] and 'S' in subMS[-1]:
            cond2_count+=1
            SEQ = SEQ[-softslen[0]:] + SEQ[:-softslen[0]]
            QUAL = QUAL[-softslen[0]:] + QUAL[:-softslen[0]]
        elif 'S' in subMS[0] and 'S' in subMS[-1]:
            cond3_count+=1
            SEQ = SEQ[-softslen[1]:] + \
                  SEQ[softslen[0]:-softslen[1]] + SEQ[:softslen[0]] 
            QUAL = QUAL[-softslen[1]:] + \
                   QUAL[softslen[0]:-softslen[1]] + QUAL[:softslen[0]]             
        #print(CIGAR, SEQ)
        outfastq.write('@'+QNAME + '\n' + SEQ + '\n+\n' + QUAL + '\n')

insam.close()
outfastq.close()

print("total count: ", total_count)
print("cond1 counts: ", cond1_count)
print("cond2 counts: ", cond2_count)
print("cond3 counts: ", cond3_count)
