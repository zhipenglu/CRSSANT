"""
gaptypes.py, Zhipeng Lu, Zhipengluchina@gmail.com

Table of Contents
0. introduction to the pipeline
1. input and output setup
2. function definitions
3. build a dictionary of non-continuous reads and output continuous ones
4. build a dictionary of segment connections and divide alignments to 4 types
5. process gapped alignments
7. process discrete chimera on the same chr and strand
8. process chimera on different chr or strand
6. process overlapping chimera
9. output counts and examples
"""







#0. introduction to the pipeline
################################################################################
################################################################################
"""
Analysis of splicing in sequencing data is very easy nowadays given the
extensive research in the field. It is more difficult to analyze gapped reads
from RNA structure data. Here we introduce a new method based on STAR, and it is
highly sensitive. This pipeline considers multiple types of alignments. This
program can process sorted or nonsorted sam files containg continuous,
gapped and chimeric alignments (see the flowchart for an outline). Gapped and
chimera are filtered based on scoreGenomicLengthLog2scale*log2(genomicLength)
and whether they match a segment dictionary. Sequencing depth affects dictionary
size, therefore a larger one can be generated from a better dataset.  

What are the major assumptions that we are making to simplify the analysis
1. Do we need to rearrange the flags? Not considered at the moment
2. multimapped chimeric alignments are not considered now (no output from STAR)
3. Only the shortest segment is checked for the penalty calculation and decision
4. Assume each alignment has at least one segment >=15nt (not necessarily true)
5. Zero-nt gaps are detectable in chimera, but not in normal arrangement.  
6. Requires large memory to store all non-continuous alignments.
7. spliced alignments, false positive short 1/2nt gaps are filtered later.
8. for discalign and diffalign not matching "connect", clips (SH) are removed 
9. calculate ligation or noncontinuous after filter: (gap1+gapm+rri+homo)/total
10. bad alignments: homopolymers, overlap with N/I, all segs < minlen


performance analysis:
1. Alignments to a dictionary: 8 sec per million, 41 mins for 300 million
2. Noncon alignments to 4 types: 112 sec per million, 73 mins for 39 million
3. Building connections: 207 mins for 498413196 connections. ~500 million
4. Gapped alignments: 132 mins for 36 million alignments.
5. Discrete chimera: 
6. Noncolinear chimera: 
7. Overlap chimera: trivial given the small number

Example discrete chimeric: NS500735:120:HH77TBGXX:1:21103:26703:6953 from file:
test1M_RfamhumanrnaMrna_gap2aChimeric.out.sam, and simplified as follows
R1 256 hs45S 10701 3 1S17M43S	 * 0 0 seq quality TAG1
R1 0   hs45S 10328 3 18S17M1N26M * 0 0 seq quality TAG2
seq:     TGGCCGTACCCATATCCGAGGGTTCCATGTGAACACAGTTGAACATGGGTCAGTCGGTCCT
quality: EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEAEEEEEEEEEEEEEEEEAE
TAG1: NH:i:2 HI:i:1 AS:i:9  nM:i:0 NM:i:0 MD:Z:17 jM:B:c,-1 jI:B:i,-1
TAG2: NH:i:2 HI:i:2 AS:i:32 nM:i:0 NM:i:0 MD:Z:43 jM:B:c,3  jI:B:i,10345,10345

basic test command:
cd /Users/lu/Documents/lulab/projects/computational/switchtest
python ~/Documents/scripts/duplex/gapanalysis.py \
test1M_RfamhumanrnaMrna_gap1bAligned.sam test1M_RfamhumanrnaMrna_gap1b -1 15

x='test30K_RfamhumanrnaMrna_gap1bAligned'
x='test1M_RfamhumanrnaMrna_gap1bAligned'
x='c'
for prefix in ${x}'cont' ${x}'gap1' ${x}'gapm' ${x}'homo' ${x}'rri' ${x}'bad';
do(samtools view -bS -o ${prefix}.bam ${prefix}.sam; \
samtools sort ${prefix}.bam ${prefix}_sorted; \
samtools index ${prefix}_sorted.bam); done

"""
################################################################################
################################################################################







#1. input and output setup
################################################################################
################################################################################
#this section sets up the input and output
#output types (5): continuous, 2-segment, multisegment, RNA-RNA inter, homorypic
#output suffix (5): cont, gap1, gapm, rri, homo
import sys, numpy, os, re, itertools, random
sys.path.append("/anaconda/lib/python2.7/site-packages")
#from intervaltree import Interval, IntervalTree
#import intervaltree
from datetime import datetime
from multiprocessing import Process, Lock, Manager

if len(sys.argv) < 6:
    print "Usage: python gapanalysis.py inputsam outprefix glenlog minlen npro"
    print "input sam can be either sorted based on coordinates or unsorted"
    print "glenlog: default -1, same as scoreGenomicLengthLog2scale"
    print "minlen: minimum length for segment to be in the database, default 15"
    print "right now the TAGs are not accurate"
    sys.exit()

inputfile = sys.argv[1]
outprefix = sys.argv[2]
glenlog = float(sys.argv[3]) #default -1. Need to test the parameters. 
minlen = int(sys.argv[4]) #min length for segment to be in the junction database
npro = int(sys.argv[5]) #number of processors to use
nonconreads = {} #dictionary to store all noncontinuous reads. need large memory
inputsam = open(inputfile, 'r'); del inputfile
inputcount = 0 #total number of reads in the file
contalign = [] #nongapped continuous alignments, due to failure in ligation
contsam = open(outprefix + 'cont.sam', 'w')
contcount = 0
gap1align = [] #all single-gap alignments (gapnorm, chimforward, chimbackward)
gap1sam = open(outprefix + 'gap1.sam', 'w') 
gap1count = 0 
gapmalign = []#all multi-gap alignments (gapnorm, chimforward, chimbackward)
gapmsam = open(outprefix + 'gapm.sam', 'w') 
gapmcount = 0
rrialign = [] #chimeric alignments on different strands or chromosomes
rrisam = open(outprefix + 'rri.sam', 'w')
rricount = 0 #some of the chimeric alignments on the same strand may be rri
homoalign = [] #negative gaps after rearrangement, likely homodimers
homosam = open(outprefix + 'homo.sam', 'w')
homocount = 0

badalign = []
badsam = open(outprefix + 'bad.sam', 'w')
badcount = 0 #bad alignments: homopolymers, chimeric with additional N/I, etc

################################################################################
################################################################################







#2. function definitions
################################################################################
################################################################################
#timenow()
#samestrand(FLAG1, FLAG2)
#getOverlap(a, b)
#mergeCIGAR(CIGAR) #get boundary of an alignment excluding S
#testconnect()
#homotypic(align1, align2) #requires re and mergeCIGAR
#filtergap1(line) #filtergap"one", filter alignment based on gap penalty
#tointerval(line) #turn an alignment to a list of intervals

#here are all the function definitions
def timenow(): return str(datetime.now())[:-7]

def samestrand(FLAG1, FLAG2):
    #check strand information based on the FLAG decimal values, e.g. 0 or 256
    if '{0:012b}'.format(int(FLAG1))[-5] == '{0:012b}'.format(int(FLAG2))[-5]:
        return True
    else: return False

def CIGARdiv(CIGAR):
    #considers all CIGAR operations [MINDSPH=X]
    #divide a cigar string to N-separated segments, used in functions below
    #return these results:
    #gaps, a list of N strings, e.g. ['100N', '200N']
    #segs, a list of segments, each with all possible operations except N
    #gaplens, a list of N lengths, e.g. [100, 200]
    #Glen, genomic length of the entire alignment (consumed ref) [MDN=X]
    #Mlens, segment lengths of matches only [M=X]
    #Rlens, segment lengths of consumed Reference [MD=X]
    #Qlens, segment lengths of consumed Query [MIS=X]
    #example: 
    #Glen,gaps,gaplens,segs,Mlens,Qlens,Rlens = CIGARdiv(CIGAR)
    #or replace unwanted output with '_' 
    gaps=re.findall('\d+N', CIGAR)
    gaplens=[int(gap[:-1]) for gap in gaps] #gap lengths 
    segs=[i.rstrip('0123456789') for i in CIGAR.split('N')]
    Mlens=[sum([int(i[:-1]) for i in re.findall('\d+[M=X]',s)]) for s in segs]
    Rlens=[sum([int(i[:-1]) for i in re.findall('\d+[MD=X]',s)]) for s in segs]
    Qlens=[sum([int(i[:-1]) for i in re.findall('\d+[MIS=X]',s)]) for s in segs]
    Glen=sum(gaplens+Rlens)
    return Glen, gaps, gaplens, segs, Mlens, Qlens, Rlens

def getOverlap(a, b):
    #determine the overlap of two intervals, both exclusive
    #returns 0 if no overlap, the overlap size if overlap
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))
    #example:print getOverlap([0,10],[10,11])

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
    #print "1S2M3N4M5I6M7S ->", mergeCIGAR("1S2M3N4M5I6M7S") #output: [1, 18]

def chimoverlap(align1, align2): #requires re and getOverlap
    #determines whether the two fragments overlap on the reference. 
    #assume the 2 alignments on the same chr and strand. Gaps included
    Rlenall1 = sum([int(i[:-1]) for i in re.findall('\d+[MDN=X]', align1[5])])
    Rlenall2 = sum([int(i[:-1]) for i in re.findall('\d+[MDN=X]', align2[5])])
    match1 = [int(align1[3]), int(align1[3]) + Rlenall1]
    match2 = [int(align2[3]), int(align2[3]) + Rlenall2]
    return getOverlap(match1, match2)

def tointerval(line): #requires CIGARdiv
    #given an alignment line, make a list of intervals, all intervals included
    align = line.split('\t')
    RNAME, POS = align[2], int(align[3])
    STRAND = '-' if '{0:012b}'.format(int(align[1]))[-5] == '1' else '+'
    _,gaps,gaplens,segs,_,_,Rlens = CIGARdiv(align[5])
    intvls = [] 
    for i in range(len(segs)): #n gaps and n+1 matches
        intvls.append([RNAME, STRAND, POS+sum(gaplens[:i])+sum(Rlens[:i]), \
                       POS+sum(gaplens[:i])+sum(Rlens[:i+1]), Rlens[i]])
    return intvls #[[RNAME, STRAND, LEFT, RIGHT, LEN], ...]


def trimseg(line, shortsegs): #requires CIGARdiv()
    #remove the short internal or external segments (do not convert to S)
    #this function does not determine whether the segment is short or not
    #the shortness of the segment is determined by the penalty and gap database
    #Here segments are defined as CIGAR substrings separated by N
    #reorganize POS, CIGAR, SEQ and QUAL only, do not change the attributes now
    align = line.split('\t')
    POS,CIGAR,SEQ,QUAL = int(align[3]),align[5],align[9],align[10]
    _,gaps,gaplens,segs,_,Qlens,Rlens = CIGARdiv(CIGAR)
    #example: gaplens [100, 200, 300, 400]
    #Qlens [5, 20, 4, 18, 6]
    #shortsegs: numbered in order? e.g. 5 segs [0,2,4]
    #keepidx, keeping index: query regions to keep, made from Qlens & shortsegs
    #Qidx: query index: [0,5,25,29,47,53]
    longsegs = [i for i in range(len(segs)) if i not in shortsegs] #e.g. [1,3]
    Qidx = [sum(Qlens[:i]) for i in range(len(segs)+1)]#e.g.[0,5,25,29,47,53]
    keepidx = [Qidx[i:i+2] for i in longsegs] #e.g. [[5,25],[29,47]]
    POSnew=POS+sum(Rlens[:longsegs[0]]+gaplens[:longsegs[0]])
    SEQnew = ''.join([SEQ[i[0]:i[1]] for i in keepidx])
    QUALnew = ''.join([QUAL[i[0]:i[1]] for i in keepidx])
    
    ##################here we remove terminal N and merge internal N for CIGAR
    gaprange = range(len(gaps)) #next make the new ops
    ops = [[segs[i],gaps[i]] if i in longsegs else gaps[i] for i in gaprange]
    if longsegs[-1]==len(gaps): ops.append(segs[-1])
    ops = sum([i if type(i) is list else [i] for i in ops], [])   
    opstart, opend = 0, len(ops)-1
    while 'N' in ops[opstart]: opstart+=1
    while 'N' in ops[opend]: opend-=1
    ops = ops[opstart:opend+1]
    opsnew = ops[:1]
    for op in ops[1:]:
        if opsnew[-1][-1]!=op[-1]: opsnew.append(op)
        else: opsnew[-1]=[str(int(opsnew[-1][:-1])+int(op[:-1]))+op[-1]]
    opsnew = sum([i if type(i) is list else [i] for i in opsnew], [])
    CIGARnew = ''.join(opsnew)
    return str(POSnew), CIGARnew, SEQnew, QUALnew

    """
    #test example: 
    linex = "A 0 chr1 50 3 5M100N20M200N4M400N18M400N6M * 0 0 \
             TGGCCGTACCCATATCCGAGGGTTCCATGTGAACACAGTTGAACATGGGTCAG \
             EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEAEEEEEEEEEE"
    shortsegs = [0,2,4]
    print trimseg(linex, shortsegs); sys.exit()
    """

def testconnect(shortsegs,longsegs,intvls,connect):
    #test which short segments cannot be matched in "connect"
    #requires the connect dictionary
    #shortsegs e.g.[0,2,4]. longsegs e.g. [1,3]
    badsegs=[]
    for i in shortsegs:
        POSshort=[k for k in range(intvls[i][2],intvls[i][3]) if not k%5]
        for j in longsegs: #only test POS that end with 0 or 5.
            POSlong=[k for k in range(intvls[j][2],intvls[j][3]) if not k%5]
            for m,n in itertools.product(POSshort, POSlong):
                if tuple(intvls[i][:2]+[m]+intvls[j][:2]+[n]) in connect: break
        else: badsegs.append(i)
    return badsegs

def trimclip(line): #this function returns the new CIGAR, SEQ and QUAL
    """requires re
    softclip fragments are a problem in rearranging CIGARs and sequences.
    Here I simply remove terminal soft/hard clips before rearranging them.
    It should only affect CIGAR, SEQ, and QUAL. 
    """
    align = line.split('\t')
    CIGAR,SEQ,QUAL = align[5],align[9],align[10]
    ops = re.findall('\d+[MINDSH=X]', CIGAR)
    CIGARnew = ''.join([i for i in ops if i[-1] not in 'SH'])
    SEQnew, QUALnew = SEQ, QUAL
    if ops[0][-1]=='S':
        SEQnew,QUALnew=SEQ[int(ops[0][:-1]):],QUAL[int(ops[0][:-1]):]
    if ops[-1][-1]=='S':
        SEQnew,QUALnew=SEQnew[:-int(ops[-1][:-1])],QUALnew[:-int(ops[-1][:-1])]
    return CIGARnew, SEQnew, QUALnew

def trimclip2line(line): #this function returns the new alignment in a string
    """ requires re
    softclip fragments are a problem in rearranging CIGARs and sequences.
    Here I simply remove terminal soft/hard clips before rearranging them.
    It should only affect CIGAR, SEQ, and QUAL. 
    """
    align = line.split('\t')
    CIGAR,SEQ,QUAL = align[5],align[9],align[10]
    ops = re.findall('\d+[MINDSH=X]', CIGAR)
    CIGARnew = ''.join([i for i in ops if i[-1] not in 'SH'])
    SEQnew, QUALnew, alignnew = SEQ, QUAL, align
    if ops[0][-1]=='S':
        SEQnew,QUALnew=SEQ[int(ops[0][:-1]):],QUAL[int(ops[0][:-1]):]
    if ops[-1][-1]=='S':
        SEQnew,QUALnew=SEQnew[:-int(ops[-1][:-1])],QUALnew[:-int(ops[-1][:-1])]
    for (i,new) in [(5,CIGARnew),(9,SEQnew),(10,QUALnew)]: alignnew[i]=new
    return '\t'.join(alignnew)

################################################################################
################################################################################







#3. build a dictionary of all reads. write out continuous alignments 
################################################################################
################################################################################
#this section processes all the header lines and builds a dictionary of all
#noncontinuous alignments (nonconreads). 
#requires a large memory to store all the non-continuous alignments
#for example, for a 100G input sam, we may need 30-50G if lots of gapped reads
#input: sam file. output: header, contsam and nonconreads dictionary

logfile = open(outprefix + 'log.out', 'w')
logstr = '\n'+timenow()+" Started gapanalysis.py ...\n"
logstr += (timenow()+" Collecting alignments into a dictionary ...\n")
logfile.write(logstr); print logstr,
header=''; firstline=''
for line in inputsam: #############process header lines for output sam files
    if line[0]=="@": header+=line
    else: firstline=line; inputcount+=1; break
contsam.write(header); gap1sam.write(header); gapmsam.write(header)
homosam.write(header); rrisam.write(header); badsam.write(header)
if "SA:Z:" not in firstline.split()[-1] and 'N' not in firstline.split()[5]:
    contalign.append(line); contcount+=1
else: nonconreads[firstline.split()[0]]=[firstline]
del firstline, header

for line in inputsam: #############store alignments in a dictionary, key=QNAME
    inputcount+=1 
    if not inputcount%10000000: #write output after every 10 million reads
        contsam.write(''.join(contalign)); contalign=[]#write cont alignments
        logstr=timenow()+" Processed "+str(inputcount)+" alignments ...\n"
        logfile.write(logstr); print logstr,
    if "SA:Z:" not in line.split()[-1] and 'N' not in line.split()[5]:
        contalign.append(line); contcount+=1; continue
    if line.split()[0] in nonconreads: nonconreads[line.split()[0]].append(line)
    else: nonconreads[line.split()[0]]=[line]
logstr=timenow()+" Total alignments: "+ str(inputcount) + '\n'
logfile.write(logstr); print logstr,
inputsam.close()
readcount=len(nonconreads)
################################################################################
################################################################################








#4. build a dictionary of all junctions and divide alignments to 4 categories
################################################################################
################################################################################
#This function divides reads (alignments) into categories and establish a
#dictionary of all good connections. return 5 dictionaries as follows. 
#this part builds a dictionary of all good gaps, similar to building a splice
#junction database, only that it is much more complicated, due to random gaps.
#input data: nonconreads, an alignment dictionary with key=QNAME
#output data include: connect,gapalign,chimdiscalign,chimoveralign,chimdiffalign

connect = {} #connect is a dictionary of all reliable connections
gapalign = {} #all normal gapped alignments
chimdiscalign = {} #chimeric alignments on same strand, chrom but no overlap
chimoveralign = {} #chimeric alignments on same strand, chrom and overlap
chimdiffalign = {} #chimeric alignments on different strand or chromosome
intvlslongall = [] #store all the long intervals to build the connections
tempcount = 0

logstr=timenow()+" Classifying alignments into 4 types ...\n"
logfile.write(logstr); print logstr,
for QNAME in nonconreads:
    tempcount+=1
    if not tempcount%1000000: #write output after every million reads
        logstr=timenow()+" Categorized "+str(tempcount)+" reads ...\n"
        logfile.write(logstr); print logstr,
    """
    ##############A. make sure that only MINSH=X are present in CIGAR
    #unnecessary to run this test for STAR output using my new parameters.
    for line in nonconreads[QNAME]: #only "MINSH=X" operations allowed
        CIGAR = line.split()[5]
        if re.findall('\d+[DP]',CIGAR): print "Unexpected:",CIGAR;sys.exit()
    """        
    ##############B. remove homopolymers artifacts. 
    SEQ = nonconreads[QNAME][0].split()[9]
    if'A'*9 in SEQ or'C'*9 in SEQ or'G'*9 in SEQ or'T'*9 in SEQ:
        badalign.append(''.join(nonconreads[QNAME])); badcount+=1; continue

    ##############C. processes gapped alignments to collect long intvls
    """
    allcigars = ''.join([i.split('\t')[5] for i in nonconreads[QNAME]])
    if "H" not in allcigars: #all of them should have 'N'
    """
    if "SA" != nonconreads[QNAME][0].split('\t')[-1][:2]: #all should have 'N'
        gapalign[QNAME] = nonconreads[QNAME]
        for line in nonconreads[QNAME]:
            intvls = tointerval(line) #get all intervals, regardless of size
            intvlslong = [intvl for intvl in intvls if intvl[4]>=minlen]
            if len(intvlslong)>=2: intvlslongall.append(intvlslong)
        continue
    
    ##############D. process chimera to collect long intvls. 
    #remaining chimera may be on diff chr or strands with more gaps
    line1, line2 = tuple(nonconreads[QNAME])
    align1, align2 = line1.split('\t'), line2.split('\t')
    if samestrand(align1[1],align2[1]) and align1[2]==align2[2]:
        if chimoverlap(align1,align2):chimoveralign[QNAME]=nonconreads[QNAME]
        else: chimdiscalign[QNAME] = nonconreads[QNAME]
    else: chimdiffalign[QNAME] = nonconreads[QNAME]
    intvls=tointerval(line1)+tointerval(line2) #[RNAME,STRAND,LEFT,RIGHT,LEN]
    intvlslong = [intvl for intvl in intvls if intvl[4]>=minlen]
    if len(intvlslong)>=2: intvlslongall.append(intvlslong)
del nonconreads


##############E. establish the connections in 5nt intervals
logstr=timenow()+" Building connections ...\n"
logfile.write(logstr); print logstr,
tempcount=0
for intvlslong in intvlslongall:
    tempcount+=1
    if not tempcount%100000:
        logstr=timenow()+" Built connections for "+str(tempcount)
        logfile.write(logstr+" alignments...\n");print logstr +" alignments..."
    grids=[] #eg.[[(chr1,'+',5),(chr1,'+',10)],[(chr2,'+',5),(chr2,'+',10)]]
    for intvl in intvlslong:
        POSs = [i for i in range(intvl[2]-4,intvl[3]+5) if not i%5]
        grids.append([(intvl[0],intvl[1],i) for i in POSs])
    for i,j in itertools.product(range(len(grids)),range(len(grids))):
        if i==j: continue
        for m,n in itertools.product(grids[i],grids[j]): connect[m+n]=''
del intvlslongall
logstr=timenow()+" Number of connections: " +str(len(connect))+'\n'
logfile.write(logstr); print logstr,

"""
#check objects at the moment
print [i for i in dir() if 'count' in i]
print [i for i in dir() if 'sam' in i]
print [i for i in dir() if 'align' in i]
print [i for i in dir() if 'temp' in i]
print [i for i in dir() if 'align' not in i and 'count' not in i \
       and 'sam' not in i and 'temp' not in i]
"""
################################################################################
################################################################################

 
 





#5. process gapped alignments
################################################################################
################################################################################
#reads with both continuous and gapped alignments are divided so only gapped are
#kept for this analysis, to speed up the algorithm.
#input data: gapalign dictionary of reads
#output files include the following: contsam, gap1sam, gapmsam

logstr=timenow()+" Processing gapped alignments ...\n"
logfile.write(logstr); print logstr,
tempcount=0

for QNAME in gapalign:
    tempcount+=1
    if not tempcount%1000000:
        logstr=timenow()+" Processed "+str(tempcount)+" gaped reads ...\n"
        logfile.write(logstr); print logstr,
    for line in gapalign[QNAME]:
        ###A. export alignments with every segment >minlen or pass penalty 
        #if there are more than one short segments, only test the shortest
        align,intvls = line.split('\t'), tointerval(line) #get all intervals
        Glen,gaps,gaplens,segs,Mlens,Qlens,Rlens = CIGARdiv(align[5])
        plusshort = glenlog*numpy.log2(Glen)+min(Mlens) #> or < 0?
        if all(map(lambda x:x[4]>=minlen,intvls)) or plusshort>=0:
            if len(intvls)==2: gap1align.append(line);gap1count+=1;continue
            if len(intvls)>2: gapmalign.append(line);gapmcount+=1;continue

        ###B. now for all remaining alignments, at least 1 seg failed the test
        shortsegs = [i for i in range(len(Qlens)) if Qlens[i]<minlen]#eg [0,2,4]
        longsegs = [i for i in range(len(Qlens)) if Qlens[i]>=minlen]#eg [1,3]
        if not longsegs:badcount+=1;badalign.append(line);continue#ignore short
        badsegs=testconnect(shortsegs,longsegs,intvls,connect)
        alignnew=align #now trim segs and output alignments
        for (i,new) in zip([3,5,9,10],trimseg(line,badsegs)):alignnew[i]=new
        linenew = '\t'.join(alignnew)
        if len(segs)-len(badsegs)<=1:contalign.append(linenew);contcount+=1
        elif len(segs)-len(badsegs)==2:gap1align.append(linenew);gap1count+=1
        else: gapmalign.append(linenew); gapmcount+=1
del gapalign
################################################################################
################################################################################







#6. process discrete chimeric alignments on the same chr and strand
################################################################################
################################################################################
#all alignments here are chimera on the same chr & strand with no overlaps.
#may be foward or backward arranged, on either strand
#input: chimdiscalign dictionary. output: contalign, gap1align, gapmalign

logstr=timenow()+" Processing discrete chimera ...\n"
logfile.write(logstr); print logstr,
chimdiscpair = {} #store the alignments to process forward and backward later
tempcount=0
for QNAME in chimdiscalign:
    tempcount+=1
    if not tempcount%1000000:
        logstr=timenow()+" Processed "+str(tempcount)+" discrete chimera ...\n"
        logfile.write(logstr); print logstr,
    
    ###A. store chimera that pass penalty test
    line1,line2 = tuple(chimdiscalign[QNAME])
    align1,align2 = line1.split('\t'),line2.split('\t')
    lineL,lineR=(line1,line2) if int(align1[3])<=int(align2[3])else(line2,line1)
    lineL,lineR = trimclip2line(lineL),trimclip2line(lineR) #remove terminal SH
    alignL,alignR = lineL.split('\t'),lineR.split('\t')
    GlenL,gapsL,gaplensL,segsL,MlensL,QlensL,RlensL = CIGARdiv(alignL[5])
    GlenR,gapsR,gaplensR,segsR,MlensR,QlensR,RlensR = CIGARdiv(alignR[5])
    Qlens = QlensL+QlensR
    gaplength = int(alignR[3])-int(alignL[3])-sum(gaplensL+RlensL)
    plusshortL = glenlog*numpy.log2(GlenL)+min(MlensL) #> or < 0?
    plusshortR = glenlog*numpy.log2(GlenR)+min(MlensR) #> or < 0?
    if len(MlensL)==1 and min(MlensL)<minlen: plusshortL=-1 #single seg < minlen
    if len(MlensR)==1 and min(MlensR)<minlen: plusshortR=-1 #single seg < minlen
    intvls = tointerval(lineL) + tointerval(lineR)
    if all(map(lambda x:x[4]>=minlen,intvls))or plusshortL>=0 and plusshortR>=0:
        chimdiscpair[QNAME] = [lineL,lineR]
        continue

    ###B. check remaining chimera against "connect", trim SH and bad segments
    #in these chimera, at least 1 segment did not pass the penalty test
    shortsegs=[i for i in range(len(Qlens)) if Qlens[i]<minlen] #eg: [0,2,4]
    longsegs=[i for i in range(len(Qlens)) if Qlens[i]>=minlen] #eg [1,3]
    badsegs = testconnect(shortsegs,longsegs,intvls,connect)
    badsegsL = [i for i in badsegs if i < len(segsL)]
    badsegsR = [i-len(segsL) for i in badsegs if i>=len(segsL)]
    goodsegsL = [i for i in range(len(segsL)) if i not in badsegsL]
    goodsegsR = [i for i in range(len(segsR)) if i not in badsegsR]

    alignnewL,alignnewR = alignL,alignR
    if goodsegsL and badsegsL:
        for (i,new) in zip([3,5,9,10],trimseg(lineL,badsegsL)):alignnewL[i]=new
    if goodsegsR and badsegsR:
        for (i,new) in zip([3,5,9,10],trimseg(lineR,badsegsR)):alignnewR[i]=new
    linenewL,linenewR = '\t'.join(alignnewL),'\t'.join(alignnewR)
    if bool(goodsegsL) != bool(goodsegsR): #one side remained, export one line
        goodsegs = goodsegsL if goodsegsL else goodsegsR
        linenew = linenewL if goodsegsL else linenewR
        if len(goodsegs)==1: contalign.append(linenew); contcount+=1
        elif len(goodsegs)==2: gap1align.append(linenew); gap1count+=1
        else: gapmalign.append(linenew); gapmcount+=1
    else: chimdiscpair[QNAME] = [linenewL,linenewR] #both sides remained
    

###now process all discrete chimera to arrange the backward or forward ones
#The bad segments have already been trimmed previously
tempcount=0
for QNAME in chimdiscpair:
    tempcount+=1
    if not tempcount%1000000:
        logstr=timenow()+" Processed "+str(tempcount)+" discrete pairs ...\n"
        logfile.write(logstr); print logstr,
            
    ##############A. set up line1/2, lineL/R, align1/2, alignL/R, CIGARL/R
    line1,line2 = tuple(chimdiscpair[QNAME])
    align1,align2 = line1.split('\t'),line2.split('\t')
    lineL,lineR=(line1,line2) if int(align1[3])<int(align2[3])else(line2,line1)
    alignL,alignR = lineL.split(),lineR.split()
    CIGARL = re.findall('\d+[MINDSH=X]', alignL[5]) #left CIGAR 
    CIGARR = re.findall('\d+[MINDSH=X]', alignR[5]) #right CIGAR

    ##############B. process foward arrangement in this step
    #only need to rearrange the CIGAR, no need to rearrange SEQ/QUAL
    #Here is an example: NS500735:120:HH77TBGXX:1:22106:7194:9151
    #hs45S 9191 4S26M33S   #hs45S 11824 30S32M1S
    revcomp = int('{0:05b}'.format(int(align1[1]))[-5])#set to 1 if rev comp
    GlenL,gapsL,gaplensL,segsL,MlensL,QlensL,RlensL = CIGARdiv(alignL[5])
    gaplength = int(alignR[3]) - int(alignL[3]) - sum(gaplensL+RlensL)
    if int(align1[3])<int(align2[3]) and not revcomp or \
       int(align1[3])>int(align2[3]) and revcomp: #forward arrangement
        CIGARF = ''.join(CIGARL+[str(gaplength)+'N']+CIGARR)
        SEQF,QUALF = alignL[9]+alignR[9],alignL[10]+alignR[10]
        alignF = align1[:5]+ [CIGARF]+align1[6:9]+[SEQF,QUALF]+align1[11:]
        gaps = re.findall('\d+N', CIGARF)
        if len(gaps)==1: gap1align.append('\t'.join(alignF)); gap1count+=1
        else: gapmalign.append('\t'.join(alignF)); gapmcount+=1
        continue

    ##############C. backward arrangement. Now CIGARs should only contain 'MIN'
    CIGARB = ''.join(CIGARL+[str(gaplength)+'N']+CIGARR)
    SEQB,QUALB = alignL[9]+alignR[9],alignL[10]+alignR[10]
    alignB = alignL[:5]+[CIGARB]+align1[6:9]+[SEQB,QUALB]+align1[11:]
    gaps = re.findall('\d+N', CIGARB)
    if len(gaps)==1: gap1align.append('\t'.join(alignB)); gap1count+=1
    else: gapmalign.append('\t'.join(alignB)); gapmcount+=1
    
del chimdiscalign, chimdiscpair
################################################################################
################################################################################






#7. process chimeric alignments on different chr or strand (non-colinear)
################################################################################
################################################################################
#this section only processes chimera on different chr or strand. min(segs)==6
#input: chimdiffalign. output: contalign, gap1align, gapmalign, rrialign
    
logstr=timenow()+" Processing noncolinear chimera (diff chr/strand) ...\n"
logfile.write(logstr); print logstr,
tempcount=0
for QNAME in chimdiffalign:
    tempcount+=1
    if not tempcount%1000000:
        logstr=timenow()+" Processed "+str(tempcount)+" noncolinear chimera ..."
        logfile.write(logstr+'\n'); print logstr
    
    ###A. all segments are at least minlen (e.g. 15nt) or pass penalty
    #only test one shortest segment, export to rri if both pass the penalty test
    line1,line2 = tuple(chimdiffalign[QNAME])
    line1,line2 = trimclip2line(line1),trimclip2line(line2) #remove terminal SH
    align1,align2 = line1.split('\t'),line2.split('\t')
    intvls = tointerval(line1)+tointerval(line2)
    Glen1,gaps1,gaplens1,segs1,Mlens1,Qlens1,Rlens1 = CIGARdiv(align1[5])
    Glen2,gaps2,gaplens2,segs2,Mlens2,Qlens2,Rlens2 = CIGARdiv(align2[5])
    Qlens = Qlens1+Qlens2
    plusshort1 = glenlog*numpy.log2(Glen1)+min(Mlens1) #>= or < 0?
    plusshort2 = glenlog*numpy.log2(Glen2)+min(Mlens2) #>= or < 0?
    if len(Mlens1)==1 and min(Mlens1)<minlen: plusshort1=-1 #single seg < minlen
    if len(Mlens2)==1 and min(Mlens2)<minlen: plusshort2=-1 #single seg < minlen
    if all(map(lambda x:x[4]>=minlen,intvls))or plusshort1>=0 and plusshort2>=0:
        rrialign.append(line1+line2);rricount+=1;continue

    ###B. now at least 1 alignment from the chimera did not pass the penalty
    #use the combined  intervals to search against the "connect" dictionary
    shortsegs=[i for i in range(len(Qlens)) if Qlens[i]<minlen] #e.g. [0,2,4]
    longsegs=[i for i in range(len(Qlens)) if Qlens[i]>=minlen] #eg [1,3]
    badsegs = testconnect(shortsegs,longsegs,intvls,connect)
    badsegs1 = [i for i in badsegs if i < len(segs1)]
    badsegs2 = [i-len(segs1) for i in badsegs if i>=len(segs1)]
    goodsegs1 = [i for i in range(len(segs1)) if i not in badsegs1]
    goodsegs2 = [i for i in range(len(segs2)) if i not in badsegs2]

    ###C. trim each alignment separately and then output together.
    alignnew1,alignnew2 = align1,align2
    if goodsegs1 and badsegs1: #segs1 are trimmed and some remained
        for (i,new) in zip([3,5,9,10],trimseg(line1,badsegs1)):alignnew1[i]=new
    if goodsegs2 and badsegs2: #segs2 are trimmed and some remained
        for (i,new) in zip([3,5,9,10],trimseg(line2,badsegs2)):alignnew2[i]=new
    linenew1,linenew2 = '\t'.join(alignnew1),'\t'.join(alignnew2)
    if bool(goodsegs1) != bool(goodsegs2): #one side remained, export one line
        goodsegs = goodsegs1 if goodsegs1 else goodsegs2
        linenew = linenew1 if goodsegs1 else linenew2
        if len(goodsegs)==1: contalign.append(linenew); contcount+=1
        elif len(goodsegs)==2: gap1align.append(linenew); gap1count+=1
        else: gapmalign.append(linenew); gapmcount+=1
    else: rrialign.append(linenew1+linenew2); rricount+=1 #both sides remained
del chimdiffalign
################################################################################
################################################################################






#8. process overlapping chimeric alignments
################################################################################
################################################################################
#input: chimoveralign. output: homoalign
#only MSH allowed in each side at the moment

logstr=timenow()+" Processing overlap chimera ...\n"
logfile.write(logstr); print logstr,
tempcount=0
for QNAME in chimoveralign:
    tempcount+=1
    if not tempcount%1000000:
        logstr=timenow()+" Processed "+str(tempcount)+" overlap chimera ...\n"
        logfile.write(logstr); print logstr,
        
    ###A. chimera on the same chr & strand, discard ones with additional 'NI'
    line1,line2 = tuple(chimoveralign[QNAME])
    align1,align2 = line1.split('\t'),line2.split('\t')
    POS1,POS2 = int(align1[3]),int(align2[3])
    lineL,lineR=(line1,line2) if POS1<=POS2 else (line2,line1)
    lineL,lineR = trimclip2line(lineL),trimclip2line(lineR) #remove terminal SH
    alignL,alignR = lineL.split('\t'),lineR.split('\t')
    POSL,POSR = int(alignL[3]),int(alignR[3])
    if 'N' in alignL[5]+alignR[5] or 'I' in alignL[5]+alignR[5]:
        badcount+=2;badalign.append(line1+line2);continue #ignore overlap w/ N/I

    ###B. convert each nt overlap to 2I1D, consuming 2nt query and 1nt reference
    Rleft =sum([int(i[:-1]) for i in re.findall('\d+[MND=X]', alignL[5])])
    Rright =sum([int(i[:-1]) for i in re.findall('\d+[MND=X]', alignR[5])])    
    overL = abs(POS1-POS2)
    overR = abs(POSR+Rright-POSL-Rleft)
    overlaplen = getOverlap([POSL,POSL+Rleft],[POSR,POSR+Rright])
    CIGARO = str(overL)+'M'+overlaplen*'2I1D'+str(overR)+'M' #overlap
    SEQO,QUALO=alignL[9]+alignR[9],alignL[10]+alignR[10]
    if POSL+Rleft>=POSR+Rright: #one side is completely within the other side
        index=[(0,overL),(overL,overL+Rright),(overL+Rright,Rleft)]
        SEQO=''.join([alignL[9][i[0]:i[1]] for i in index[0:2]+index[1:3]])
        QUAL=''.join([alignL[10][i[0]:i[1]] for i in index[0:2]+index[1:3]])
    for (i,O) in [(5,CIGARO),(9,SEQO),(10,QUALO)]: alignL[i]=O #modify alignL
    homoalign.append('\t'.join(alignL)); homocount+=1
del chimoveralign
#Note: chimoveralign: intermediate dictionary for potential chimeric overlapped
#homoalign: final chimeric overlapped alignments (homotypic) to be output
################################################################################
################################################################################







#9. output and examples
################################################################################
################################################################################
#this section is output and examples
for x,y in[(contalign,contsam),(gap1align,gap1sam),(gapmalign,gapmsam),\
           (rrialign,rrisam),(homoalign,homosam),(badalign,badsam)]:
    y.write(''.join(x)); del x; y.close() #write to output
logstr=timenow()+" Finished gapanalysis.py successfully\n\n"
logfile.write(logstr); print logstr,

logstr= \
"              Total input alignment number: " + str(inputcount) + '\n' + \
"           Continuous alignments (no gaps): " + str(contcount) + '\n' + \
"             Two-segment gapped alignments: " + str(gap1count) + '\n' + \
"           Multi-segment gapped alignments: " + str(gapmcount) + '\n' + \
"        Other chimeric (different str/chr): " + str(rricount) + '\n' + \
"          Overlapping chimeric (homotypic): " + str(homocount) + '\n' + \
"                            Bad alignments: " + str(badcount) + '\n'
logfile.write(logstr); print logstr
logfile.close()



""" Example output:
              Total input alignment number: 278030
           Continuous alignments (no gaps): 199953
             Two-segment gapped alignments: 42442
           Multi-segment gapped alignments: 559
        Other chimeric (different str/chr): 11095
          Overlapping chimeric (homotypic): 118
                            Bad alignments: 767

continuous alignments from the raw file: 
from gapalign:      cont: 5701  gap1:31083  gapm:555
from chimdiscalign: cont: 0     gap1:11359  gapm:4    
from chimdiffalign: cont: 1                          rri:11095
from chimoveralign:                                  homo: 118
"""

"""
#Using the following script check output from each processing section:
print "input", inputcount
print " cont", contcount
print " gap1", gap1count
print " gapm", gapmcount
print "  rri", rricount
print " homo", homocount
print "  bad", badcount
for x,y in[(contalign,contsam),(gap1align,gap1sam),(gapmalign,gapmsam),\
           (rrialign,rrisam),(homoalign,homosam),(badalign,badsam)]:
    y.write(''.join(x)); del x; y.close() #write to output
print timenow(),"So far so good."; sys.exit()
"""
################################################################################
################################################################################
