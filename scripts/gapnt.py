"""
gapnt.py. Zhipeng Lu, zhipengluchina@gmail.com, 2019-11-30
This script analyzes the gap nt identity of the gapped reads data from other
methods that use psoralen crosslinking to determine RNA duplexes.

This script is relatively fast on the Rfam+mRNA collections, but super slow on
the whole genomes due to the larger genome and insertions.


test command
cd /Users/lu/Documents/lulab/projects/computational/othermethods

LIGR example. To analyze individual RNAs
awk '$3=="ACTB"' SRR3361013_RfamhumanrnaMrna_gap1Aligned_N_prim.sam > \
SRR3361013_RfamhumanrnaMrna_gap1Aligned_prim_N_ACTB.sam
python ~/Documents/scripts/duplex/gapnt.py ref_ACTB.fa \
SRR3361013_RfamhumanrnaMrna_gap1Aligned_prim_N_ACTB.sam

LIGR example. To analyze all reads: 
python ~/Documents/scripts/duplex/gapnt.py \
~/Documents/chang/psoralen/interaction/Rfam/RfamhumanrnaMrna.fa \
SRR3361013_RfamhumanrnaMrna_gap1Aligned_prim_N.sam

LIGR example. To analyze all RNA except hs45S and the Rfam collection 
awk '($3!="hs45S")&&($3!~/^RF0/)' \
SRR3361013_RfamhumanrnaMrna_gap1Aligned_N_prim.sam > \
SRR3361013_RfamhumanrnaMrna_gap1Aligned_N_prim_mRNA.sam #mostly mRNAs
python ~/Documents/scripts/duplex/gapnt.py \
~/Documents/chang/psoralen/interaction/Rfam/RfamhumanrnaMrna.fa \
SRR3361013_RfamhumanrnaMrna_gap1Aligned_prim_N_mRNA.sam

LIGR example. To visualize all the alignments with 1 or 2nt gaps: 
awk '($1~/^@/)||($6~/M1N/)||($6~/M2N/)' \
SRR3361013_RfamhumanrnaMrna_gap1Aligned_N_prim.sam > \
SRR3361013_RfamhumanrnaMrna_gap1Aligned_1N2N_prim.sam

Analysis showed a dramatic bias towards T, consistent with psoralen reaction
It is more obvious for the mRNAs, where other types of mistakes are minimal.
Noncoding RNAs have lots of paralogs and repeats that could produce errors. 
ACTB TUBB RPL10 RPL_RPS RF00100_7SK mRNA_1or2nt  mRNA_1nt
A 6  13   39    4995    5722        36532  19069        
C 44 52   28    4185    10259       29996  14266        
G 16 54   88    5145    6210        40317  18309        
T 53 79   124   12965   9060        86594  43429        


#######################The following is analysis of SPLASH data
python ~/storage/zhipengl/bin/gapnt.py \
~/storage/zhipengl/staridx/starRfamhumanrnaMrna/RfamhumanrnaMrna.fa \
SRR3404937_RfamhumanrnaMrna_gap1Aligned_prim_N.sam

#######################The following is analysis of PARIS data
python ~/storage/zhipengl/bin/gapnt.py \
~/storage/zhipengl/staridx/starRfamhumanrnaMrna/RfamhumanrnaMrna.fa \
AMT_Stress_trim_nodup_bc07_RfamhumanrnaMrna_gap1Aligned_prim_N.sam
"""

import sys, re
from datetime import datetime

if len(sys.argv)<3:
    print "\nUsage: python gapnt.py fasta inputsam"
    print "fasta is the reference sequence"
    print "distribution of nt is printed to the stdout\n"
    sys.exit()

fasta = open(sys.argv[1], 'r')
inputsam = open(sys.argv[2], 'r')
gap1s = '' #string to store all nts from 1nt gaps
gap2s = '' #string to store all nts from 2nt gaps
gap3s = '' #string to store all nts from 2nt gaps
gapms = '' #string to store all nts from multiple nt (>2) gaps

#build a dictionary of all fasta references.
fastadic={}
linecount=0
for line in fasta:
    linecount+=1
    if not linecount%1000000:
        print str(datetime.now())[:-7], "Processed",linecount,"lines fasta ..."
    if line[0]=='>': RNAME=line.split()[0].strip('>\n'); fastadic[RNAME]=[]
    else: fastadic[RNAME].append(line.strip('\n'))
for RNAME in fastadic: fastadic[RNAME]=''.join(fastadic[RNAME])
print "Number of references:", len(fastadic)

        

#for each line, examine the gap length and extract nt frequency
tempcount=0
for line in inputsam:
    if line[0] == '@': continue
    tempcount+=1
    if not tempcount%1000000:
        print str(datetime.now())[:-7], "Processed",tempcount,"lines sam ..."
    align = line.split()
    RNAME, POS, CIGAR = align[2], int(align[3]), align[5]
    refcon = re.findall('\d+[MDN=X]', CIGAR) #CIGAR operations that consume ref
    #refcon: 3S10M4N20M1I6M1N9M convert to ['10M','4N','20M','6M','1N','9M']
    gaplens = [int(i[:-1]) for i in refcon if 'N' in i]
    conlist = [['', 0, 0]]
    for i in refcon:conlist+=[[i[-1],conlist[-1][2],conlist[-1][2]+int(i[:-1])]]
    #conlist: ['10M','4N','20M','6M','1N','9M'] to [['M',0,10],['N',10,14], ...]
    gap1loc=[[i[1],i[2]] for i in conlist[1:] if 'N'==i[0] and i[2]-i[1]==1]
    gap2loc=[[i[1],i[2]] for i in conlist[1:] if 'N'==i[0] and i[2]-i[1]==2]
    gap3loc=[[i[1],i[2]] for i in conlist[1:] if 'N'==i[0] and i[2]-i[1]==3]
    gapmloc=[[i[1],i[2]] for i in conlist[1:] if 'N'==i[0] and i[2]-i[1]>3 and \
             i[2]-i[1]<11]

    for i in gap1loc: gap1s+=fastadic[RNAME][POS+i[0]-1:POS+i[1]-1].upper()
    for i in gap2loc: gap2s+=fastadic[RNAME][POS+i[0]-1:POS+i[1]-1].upper()
    for i in gap3loc: gap3s+=fastadic[RNAME][POS+i[0]-1:POS+i[1]-1].upper()
    for i in gapmloc: gapms+=fastadic[RNAME][POS+i[0]-1:POS+i[1]-1].upper()

ntfreq1=[gap1s.count('A'),gap1s.count('C'),gap1s.count('G'), gap1s.count('T')]
ntfreq2=[gap2s.count('A'),gap2s.count('C'),gap2s.count('G'), gap2s.count('T')]
ntfreq3=[gap3s.count('A'),gap3s.count('C'),gap3s.count('G'), gap3s.count('T')]
ntfreqm=[gapms.count('A'),gapms.count('C'),gapms.count('G'), gapms.count('T')]

print "Frequency of nts for 1nt, 2nt, 3nt and multi-nt gaps, in the order ACGT:"
print ntfreq1, '\n', ntfreq2, '\n', ntfreq3, '\n', ntfreqm, '\n'
inputsam.close()
print str(datetime.now())[:-7],"Finished gapnt.py analysis successfully.\n"


