#!/usr/bin/python3
import sys, numpy, os, re, itertools, random
import pysam

sys.path.append("/anaconda/lib/python2.7/site-packages")
#from intervaltree import Interval, IntervalTree
#import intervaltree
from datetime import datetime
from multiprocessing import Process, Lock, Manager

if len(sys.argv) < 6:
    print("Usage: python gapanalysis.py inputbam outputbam glenlog minlen npro")
    print("input sam can be either sorted based on coordinates or unsorted")
    print("glenlog: default -1, same as scoreGenomicLengthLog2scale")
    print("minlen: minimum length for segment to be in the database, default 15")
    print("right now the TAGs are not accurate")
    sys.exit()

inputfile = sys.argv[1]
outprefix = sys.argv[2]
glenlog = float(sys.argv[3]) #default -1. Need to test the parameters. 
minlen = int(sys.argv[4]) #min length for segment to be in the junction database
npro = int(sys.argv[5]) #number of processors to use
nonconreads = {} #dictionary to store all noncontinuous reads. need large memory

starttime = datetime.now()
if(inputfile.endswith("sam")):
    print(inputfile," is SAM")
    inputbam = pysam.AlignmentFile(inputfile, "r", require_index=False)
if(inputfile.endswith(".bam")):
    print(inputfile, " is BAM")
    inputbam = pysam.AlignmentFile(inputfile, "rb", require_index=False)

if(inputbam is None):
    print("Usage: gaptypes3.py inputbam outprefix glenlog minlen npro")
    print("   input: sam/bam file can be either sorted based on coordinates or unsorted. must use .sam or .bam input file")
    print("   glenlog: default -1, same as scoreGenomicLengthLog2scale")
    print("   minlen: minimum length for segment to be in the database, default 15")
    print("   right now the TAGs are not accurate")

    sys.exit()


inputcount = 0 #total number of reads in the file
readcount = 0
continuous_count = 0
noncontinuous_count = 0;

outputbam = pysam.AlignmentFile(outprefix+".bam", "wb", template=inputbam)
logfile = open(outprefix + 'log.out', 'w')

max_reads = 3000000
max_aligns = 10000000


################################################################################

def timenow(): return str(datetime.now())[:-7]

def logmsg(logstr):
    logfile.write(logstr); 
    print(logstr, end='')




################################################################################
#
# main loop write out continuous alignments and group reads for analysis
#
################################################################################

#might need to do 2 passes through the file
# 1st pass to find the intvlslongall and connections hashes (should be small memory)
# 2nd pass to process read QNAME alignment-groups


#pass #1 performs check for continuous vs non-continuous and build long-connections
inputcount=0
readcount=0
badcount=0
current_qname = ""
current_aligns = []
t_noncont_count=0;
for align in inputbam:
    inputcount+=1 
    
#    if not inputcount%1000000: #write output after every 1 million reads
#        t_time = datetime.now()
#        difftime = (t_time-starttime).seconds
#        rate = ((inputcount/1000000) / difftime) * 60
#        logmsg(timenow()+" Phase1-connections "+str(inputcount/1000000)+"mil alignments ["+
#               str(readcount)+" reads] "+
#               "("+str(difftime)+"sec, "+("%.3f" % rate)+" mega-align/min)\n")

    QNAME = align.query_name  #line[0]
    CIGAR = align.cigarstring  #line[5]
    FLAG = align.flag
    #SEQ   = align.query_sequence
    #QUAL  = align.qual

    # this is the only logic for deciding if line goes into nonconreads
    #if "SA:Z:" not in line.split()[-1] and 'N' not in line.split()[5]:
    #    contalign.append(line); contcount+=1; continue
    #if line.split()[0] in nonconreads: nonconreads[line.split()[0]].append(line); noncontinuous_count+=1
    #else: nonconreads[line.split()[0]]=[line]; noncontinuous_count+=1


    if((not align.has_tag('SA')) and ('N' not in CIGAR)):
        continuous_count+=1; 
    else:
        noncontinuous_count+=1
    
    if(current_qname==""): 
        #print("first qname ",QNAME)
        current_qname = QNAME
        current_aligns = []
        readcount+=1

    if(QNAME != current_qname):
        for t_align in current_aligns: 
            outputbam.write(t_align);            
        current_qname = QNAME  #reset to the next QNAME block
        if(readcount >= max_reads): break
        #if(inputcount >= max_aligns): break
        del current_aligns
        
        if(readcount%100000==0): #write output after every 1 million reads
            t_time = datetime.now()
            difftime = (t_time-starttime).seconds
            rate = (inputcount / difftime)
            logmsg("read input "+ f"{inputcount:,}" +" alignments ["+ f"{readcount:,}" +" reads] "+
                  "("+str(difftime)+"sec, "+ f"{rate:,}" +" align/sec)\n")

        current_aligns=[] #release all the align from current_qname to save memory
        readcount+=1

    current_aligns.append(align)


endtime = datetime.now()
#print("finished reading inputbam in %1.3f seconds\n" % ((endtime-starttime).seconds))
logstr= str("finished read inputbam in "+ ("%1.3f"%((endtime-starttime).seconds)) +" seconds\n" + 
            "                         input reads: " + f"{readcount:,}" +"\n" + 
            "              input alignment number: " + f"{inputcount:,}" + '\n' + 
            "     Continuous alignments (no gaps): " + f"{continuous_count:,}" + '\n' + 
            "           Non-Continuous alignments: " + f"{noncontinuous_count:,}" + '\n')

logmsg(logstr)

inputbam.close()
outputbam.close()
logfile.close()



