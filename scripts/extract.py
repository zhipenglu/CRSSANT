#!/usr/bin/python3
import sys, numpy, os, re, itertools, random
import pysam
import argparse

sys.path.append("/anaconda/lib/python2.7/site-packages")
#from intervaltree import Interval, IntervalTree
#import intervaltree
from datetime import datetime
from multiprocessing import Process, Lock, Manager

max_reads = 3000000
max_aligns = 10000000
mode = "maxreads"


parser = argparse.ArgumentParser(
                    prog='extract.py',
                    description='extract alignments from a QNAME (-n) sorted bam/sam based on either maxread or maxalignment counts. Always returns all aligments for each read')

parser.add_argument('inputfile', help="sam/bam filename")           # positional argument
parser.add_argument('outprefix', help="output prefix" )         # positional argument
parser.add_argument('-r', '--maxreads', type=int,
                    help="stop extraction at max read count <default 3000000 reads>")
parser.add_argument('-a', '--maxaligns', type=int,
                    help="stop extraction when alignment count exceeds max. always extracts all alignments for each read")
parser.add_argument('-v', '--verbose', action='store_true')  # on/off flag

args = parser.parse_args()
#print(args)

inputfile = args.inputfile
outprefix = args.outprefix

if(args.maxreads != None):
    max_reads = args.maxreads
    mode = "maxreads"
    print("mode maxreads",max_reads)
elif(args.maxaligns != None):
    max_aligns = args.maxaligns
    mode = "maxaligns"
    print("mode alignments",max_aligns)


starttime = datetime.now()
inputbam = None
if(inputfile.endswith("sam")):
    print(inputfile," is SAM")
    try:
        inputbam = pysam.AlignmentFile(inputfile, "r", require_index=False)
    except OSError as e:
        sys.exit()
if(inputfile.endswith(".bam")):
    print(inputfile, " is BAM")
    try:
        inputbam = pysam.AlignmentFile(inputfile, "rb", require_index=False)
    except OSError as e:
        sys.exit()
if(inputbam is None): 
    print("\nERROR: input file incorrect\n")
    parser.print_usage() #print_help
    sys.exit()


inputcount = 0 #total number of reads in the file
readcount = 0
continuous_count = 0
noncontinuous_count = 0;

outputbam = pysam.AlignmentFile(outprefix+".bam", "wb", template=inputbam)
logfile = open(outprefix + 'log.out', 'w')


################################################################################

def timenow(): return str(datetime.now())[:-7]

def logmsg(logstr):
    logfile.write(logstr); 
    print(logstr, end='')


################################################################################
#
# main loop group alignments by read(QNAME) and write out until max cutoff
#
################################################################################

inputcount=0
readcount=0
badcount=0
current_qname = ""
current_aligns = []
t_noncont_count=0;
for align in inputbam:
    inputcount+=1 
    
    QNAME = align.query_name  #line[0]
    CIGAR = align.cigarstring  #line[5]
    FLAG = align.flag
    #SEQ   = align.query_sequence
    #QUAL  = align.qual

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

        if(mode=="maxreads" and (readcount >= max_reads)): break
        if(mode=="maxaligns" and (inputcount >= max_aligns)): break

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



