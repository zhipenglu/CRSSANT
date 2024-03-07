#!/usr/bin/python3

"""
gapfilter.py, Zhipeng Lu. 2020, zhipengluchina@gmail.com. 
This script filters PARIS or other data to remove splice junction alignments
and short deletions, e.g. 1-2nt gaps.

Get junctions from a GTF file and compare to non-continuous alignments.
All CIGAR operations are considered to ensure compatibility with other input SAM
The splicing junctions were shifted to the left and right by 2 nt to account for
gap openning errors. 
"""
import sys, os, re
from datetime import datetime
import pysam
import argparse

parser = argparse.ArgumentParser(
                    prog='gapfilter.py',
                    description='Output files gap1.bam and gapm.bam may contain alignments that have only splicing junctions and short 1-2 nt gaps due to artifacts. These are filtered out using gapfilter.py before further processing.')

parser.add_argument('annofile', help="GTF annotation filename")  # positional argument
parser.add_argument('inputfile', help="sam/bam input filename")  # positional argument
parser.add_argument('outprefix', help="output prefix" )          # positional argument
parser.add_argument('--idloc', type=int, default=11,
                    help="space sep. column num of transcript_id in GTF, (default 11)")
parser.add_argument('--short', action='store_true',
                    help="remove short 1-2nt gaps (default is to 'ignore short 1-2nt gaps')")
parser.add_argument('-v', '--verbose', action='store_true')  # on/off flag

args = parser.parse_args()

annofile = args.annofile
inputfile = args.inputfile
outprefix = args.outprefix
idloc = args.idloc
short = args.short #boolean: true mean remove alignments with only 1/2nt gaps

starttime = datetime.now()

logfile = open(outprefix + '_log.out', 'w')
def logmsg(*argv):
  logstr = " ".join(str(e) for e in argv)
  logfile.write(logstr+"\n"); 
  print(logstr)

logmsg("annofile:  ", annofile)
logmsg("inputfile: ", inputfile)
logmsg("outprefix: ", outprefix)
logmsg("idloc:     ", idloc)
logmsg("short:     ", short)

#if len(sys.argv) < 6:
#    print("\nUsage: python gapfilter.py annotation insam outsam idloc short\n"
#     "annotation should be in the GTF format\n"
#     "Find location of the transcript_id (space sep. column num) for idloc\n"
#     "idloc=11 for hg38_refGene.sorted.gtf\n"
#     "shift of 2nt is used to account for mismapping to splice junctions\n"
#     "optional: remove short (1/2nt) gaps, yes or no\n"
#     "the gapmfilter.sam output contain alignments with 1 or more good gaps\n")
#    sys.exit()

anno = open(annofile, 'r')

inputbam = None
if(inputfile.endswith("sam")):
    logmsg(inputfile," is SAM")
    try:
        inputbam = pysam.AlignmentFile(inputfile, "r", require_index=False)
    except OSError as e:
        sys.exit()
if(inputfile.endswith(".bam")):
    logmsg(inputfile," is BAM")
    try:
        inputbam = pysam.AlignmentFile(inputfile, "rb", require_index=False)
    except OSError as e:
        sys.exit()
if(inputbam is None): 
    logmsg("\nERROR: input file incorrect\n")
    parser.print_usage() #print_help
    sys.exit()

outputbam = pysam.AlignmentFile(outprefix+".bam", "wb", template=inputbam)

logmsg(str(datetime.now())[:-7], "Starting the gapfilter analysis ...")

shift = 2 #allow the junction to shift to the left or right by at most 2nt. 
inputcount = 0
gapcount = 0
ligcount = 0
ligmcount =0 # number of alignments with more than 1 good gaps. 




def getgaps(align): #turn an alignment into gaps [(RNAME, STRAND L, R), ...]
    #align = line.split()
    #RNAME, POS, CIGAR = align[2], int(align[3]), align[5]
    #STRAND = '-' if '{0:012b}'.format(int(align[1]))[-5] == '1' else '+'
    RNAME  = inputbam.get_reference_name(align.reference_id)
    POS    = align.reference_start +1  #line[3] GTF and SAM are 1 based
    CIGAR  = align.cigarstring  #line[5]
    STRAND = '-' if (align.is_reverse) else '+'

    gaps = [] #store all gaps from this CIGAR string, each as a 3-tuple. 
    Ns = [int(i[:-1]) for i in re.findall('[0-9]+N', CIGAR)] #gap lengths
    arms =[i.rstrip('0123456789') for i in CIGAR.split('N')]
    rx = [] #reference consumed: MD=X
    for arm in arms:
        rx.append(sum([int(i[:-1]) for i in re.findall('[0-9]+[MD=X]', arm)]))
    for i in range(len(Ns)): #combine ref and gap lengths to make the junctions
        l, r = POS+sum(rx[:i+1])+sum(Ns[:i])-1, POS+sum(rx[:i+1])+sum(Ns[:i+1])
        gaps.append((RNAME, STRAND, l, r))
    return gaps
    

#############Process all GTF data into a dictionary of transcripts
#GTF 9 fields: seqname source feature start end score strand frame attribute
logmsg(str(datetime.now())[:-7], "Reading the annotation GTF file ...")
transdict = {} #transcript dictionary, key: transcript id, value: exon bounds 
prev_transcript = ""
for line in anno:
    record = line.split()
    if line[0]=="#" or record[0]=='track' or record[2]!='exon': continue
    if record[idloc-1] != "transcript_id": #9-12: gene_id "X";transcript_id "Y";
        logmsg("transcript_id not in expected location, exiting"); sys.exit()
    transcript = record[idloc].strip('";')
    if(prev_transcript==""): prev_transcript=transcript
    coord = [record[0], record[6], int(record[3]), int(record[4])]
    if transcript not in transdict: transdict[transcript] = [coord]
    else: transdict[transcript].append(coord)
    #if prev_transcript != transcript: logmsg(prev_transcript, transdict[prev_transcript])
    prev_transcript = transcript
#print(transdict)


#############Store all the junctions in a dictionary, each junction as a tuple
#allow shifts to the left or right, e.g. by 2nt (set by the shift option). 
logmsg(str(datetime.now())[:-7], "Building the splicing junction database ...")
junctdict = {} #key: (seqname, exon1_end, exon2_start)
for transcript in transdict:
    transdict[transcript]=sorted(transdict[transcript])
    exons = len(transdict[transcript])
    if exons == 1: continue
    for i in range(exons-1):
        exon = transdict[transcript][i]
        for j in range(-shift, shift+1): 
            junctdict[(exon[0], exon[1], exon[3]+j, \
                       transdict[transcript][i+1][2]+j)] =''


logmsg(str(datetime.now())[:-7], "Started filtering splice junctions ...")
#############Find locations of all gaps, only N
#tolerate ambiguities at the junctions if off by 2nt as defined by shift

for align in inputbam:
    CIGAR = align.cigarstring  #line[5]

    inputcount +=1
    if not inputcount%1000000:
        logmsg(str(datetime.now())[:-7], "Processed", inputcount, "reads ...")
    
    if 'N' not in CIGAR: continue
    gaps = getgaps(align); gapcount +=1

    #short gaps are tolerated.
    goodidex = [0]
    if not short: #only remove alignments that only have splice junctions.
        goodidx=[0 if gap in junctdict else 1 for gap in gaps]
    elif short: #remove alignments where all gaps are short or spliced
        goodidx=[0 if gap in junctdict or gap[3]-gap[2]<4 else 1 for gap in gaps]
    if max(goodidx)==1: #0: bad gaps (spliced or 1-2nt); 1: good gaps 
        #at least one gap is good. 
        outputbam.write(align);
        ligcount+=1 
        if len([i for i in goodidx if i==1])>1: ligmcount+=1
        

logmsg("\n                        Total alignments:", inputcount)
logmsg("                   All gapped alignments:", gapcount)
logmsg("    Alignments with at least 1 good gaps:", ligcount)
logmsg("    Alignments with at least 2 good gaps:", ligmcount)
logmsg("  Number of annotated splicing junctions:", len(junctdict)/5, '\n')
logmsg(str(datetime.now())[:-7],"Finished gapfilter.py successfully")

anno.close()
inputbam.close()
outputbam.close()

runtime = (datetime.now()-starttime).seconds
if(runtime>3600):
  logmsg("total runtime "+ f"{(runtime/3600.0):,.1f}" + " hours")  
elif(runtime <200):
  logmsg("total runtime "+ f"{runtime:,.1f}" + " seconds")
else:
  logmsg("total runtime "+ f"{(runtime/60.0):,.1f}" + " minutes")
logfile.close()

"""
Analysis of AMT_Stress_trim_nodup_bc07_starhg38Aligned_prim_N.nosplice.sam
                        Total alignments: 3570138
                   All gapped alignments: 3570138
Alignments with at least 1 nonsplice gap: 3549243
  Number of annotated splicing junctions: 261919
It seems that spliced reads are not a major problem in my final output. 
"""
