"""
dgsim.py

Generated random numbers, and used initialized seeds to maintain consistency.

This script generates a simple set of alignments in DGs to test crssant.py
reference: chr1:1-1000 (first 10kb of hg38 chr1 happens to be a stretch of N).
DG core coordinates range from 100 to 900 (leave out 100 from both ends)
one gene: GENE1 end to end
No. of DGs: 100. Alignments in each DG: 10 to 1000 (variable).
Each arm is defined by a core and extends to each side (both variable).
No two DGs should overlap on the core regions of both arms.

Goals: 
1. To test processing speed: vary reads numbers. 
2. To test read length: vary edge length
3. To test overlap threshold t_o: vary core length
4. To test the two clustering method. 

Files: 
chr1.fa
chr1_len.txt
chr1_gene.bed: chr1 0 1000 GENE1 1000 +

Steps:
1. Generate a list of DG cores (format: 4-tuple), with a predefined core length
2. Convert the DG core list to a dictionary of DGs by extending on each side
3. Convert the DGs with intervals to a sam file. 
4. Test crssant.py on the sam file. 



Example commands:
cd ~/Documents/scripts/crssant
python dgsim.py chr1_gap1sim_core10.sam
samtools view -bS -o chr1_gap1sim_core10.bam chr1_gap1sim_core10.sam
samtools sort chr1_gap1sim_core10.bam chr1_gap1sim_core10_sorted
samtools index chr1_gap1sim_core10_sorted.bam
bedtools genomecov -bg -split -strand + -ibam chr1_gap1sim_core10_sorted.bam -g \
chr1_len.txt > chr1_gap1sim_core10_plus.bedgraph
bedtools genomecov -bg -split -strand - -ibam chr1_gap1sim_core10_sorted.bam -g \
chr1_len.txt > chr1_gap1sim_core10_minus.bedgraph
"""




import sys
import random
import itertools
if len(sys.argv) < 2:
    print("Usage: python dgsim.py sim.sam")
    print("Parameters set within this script. Output is sim.sam")



################################################################################
#####custom parameters, can be varied to test crssant performance. 
RNAME = 'chr1' #reference name
RLEN = 500 #reference genome length
edge = 100 #the 5' and 3' ends excluded from the DGcore random sampling
corelen = 5 #test 1 to 10.
exlower, exupper = 5, 15 #lower and upper bounds of extra length on each side
coregap = 50 #gap length between the two core regions
corewithin = 100 #
numDGs = 100 #number of simulated DGs
DGlower, DGupper = 10,100 #lower and upper bounds of read number in each DG
header = '@HD\tVN:1.4\n@SQ\tSN:' + RNAME + '\tLN:' + str(RLEN) + '\n'
TAGs = 'NH:i:1\tHI:i:1\tAS:i:0\tnM:i:0\tNM:i:0\tMD:Z:0\tjM:B:c,0\tjI:B:i,0,0'



################################################################################
#####generate a list of DG cores
DGcores = [] #list of DG cores
DGcoresxy = {} #2-tuples specifying all possible connections between the 2 arms
seed = 0 #random seed
while len(DGcores) < numDGs:
    seed+=1
    random.seed(seed)    
    i=random.randrange(edge,RLEN-edge)
    random.seed(seed+numDGs)
    j=random.randrange(edge,RLEN-edge)
    if coregap<j-i<corewithin:
        DGcore = [i,i+corelen,j,j+corelen]
        overlap = 0
        combi = itertools.product(range(i,i+corelen), range(j,j+corelen))
        for (x,y) in combi:
            if (x,y) in DGcoresxy: overlap += 1
        if overlap == 0:
            for (x,y) in combi: DGcoresxy[(x,y)] = 0
            DGcores.append(DGcore)
#print("DGcores:", DGcores)
################################################################################




################################################################################
#####produce a dictionary of DGs, key = DGcore (4-tuple), value = DGintvls
DGdict = {} 
numaligns = [] #numbers of alignments in each DG           
for i in range(numDGs):
    random.seed(i)
    numalign = random.randint(DGlower,DGupper)
    numaligns.append(numalign)
    DGintvls = []
    DGcore = DGcores[i]
    for x in range(numalign):
        random.seed(x*numalign+1); rnd0 = random.randint(exlower,exupper)
        random.seed(x*numalign+2); rnd1 = random.randint(exlower,exupper)
        random.seed(x*numalign+3); rnd2 = random.randint(exlower,exupper)
        random.seed(x*numalign+4); rnd3 = random.randint(exlower,exupper)
        DGintvl = [DGcore[0]-rnd0,DGcore[1]+rnd1,DGcore[2]-rnd2,DGcore[3]+rnd3]
        DGintvls.append(DGintvl)
    DGdict[tuple(DGcore)] = DGintvls
#print("Numbers of alignments in each DG:", numaligns)
#print(DGdict.keys())
################################################################################




################################################################################
#####convert the DG dictionary into a sam file of alignments
#QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL TAGs
outsam = open(sys.argv[1], 'w')
outsam.write(header)
for DG in DGdict:
    dgid = 'core_' + '_'.join([str(i) for i in DG])
    #print(dgid)
    for i in range(len(DGdict[DG])):
        intvl = DGdict[DG][i]
        QNAME = dgid + '_read_' + str(i)
        FLAG = '0'
        POS = str(intvl[0])
        MAPQ = '255'
        M1,N1,M2 = intvl[1]-intvl[0],intvl[2]-intvl[1],intvl[3]-intvl[2]
        CIGAR = str(M1)+'M'+str(N1)+'N'+str(M2)+'M'
        SEQ = (M1+M2)*'N'
        QUAL = (M1+M2)*'A'
        align = [QNAME,FLAG,RNAME,POS,MAPQ,CIGAR,'*','0','0',SEQ,QUAL,TAGs]
        line = '\t'.join(align) + '\n'
        outsam.write(line)
outsam.close()
################################################################################





