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
parser.add_argument('-s', '--short', action='store_true',
                    help="remove short 1-2nt gaps (default is to 'ignore short 1-2nt gaps')")
parser.add_argument('-v', '--verbose', action='store_true')  # on/off flag

args = parser.parse_args()

annofile = args.annofile
inputfile = args.inputfile
outprefix = args.outprefix
#inputfile = sys.argv[2]
#outprefix = sys.argv[3]
idloc = args.idloc
short = args.short #boolean: true mean remove alignments with only 1/2nt gaps

print("annofile:  ", annofile)
print("inputfile: ", inputfile)
print("outprefix: ", outprefix)
print("idloc:     ", idloc)
print("short:     ", short)


#if len(sys.argv) < 6:
#    print("\nUsage: python gapfilter.py annotation insam outsam idloc short\n"
#     "annotation should be in the GTF format\n"
#     "Find location of the transcript_id (space sep. column num) for idloc\n"
#     "idloc=11 for hg38_refGene.sorted.gtf\n"
#     "shift of 2nt is used to account for mismapping to splice junctions\n"
#     "optional: remove short (1/2nt) gaps, yes or no\n"
#     "the gapmfilter.sam output contain alignments with 1 or more good gaps\n")
#    sys.exit()

#anno = open(sys.argv[1], 'r')
#inputsam = open(sys.argv[2], 'r')
#outputsam = open(sys.argv[3], 'w')
#idloc = int(sys.argv[4])
#short = sys.argv[5] #yes to remove alignments with only 1/2nt gaps, or no. 

anno = open(annofile, 'r')

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

outputbam = pysam.AlignmentFile(outprefix+".bam", "wb", template=inputbam)
#logfile = open(outprefix + 'log.out', 'w')

print(str(datetime.now())[:-7], "Starting the gapfilter analysis ...")

shift = 2 #allow the junction to shift to the left or right by at most 2nt. 
inputcount = 0
gapcount = 0
ligcount = 0
ligmcount =0 # number of alignments with more than 1 good gaps. 
outstring = ''


def getgaps(line): #turn an alignment into gaps [(RNAME, STRAND L, R), ...]
    align = line.split()
    RNAME, POS, CIGAR = align[2], int(align[3]), align[5]
    STRAND = '-' if '{0:012b}'.format(int(align[1]))[-5] == '1' else '+'
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
print(str(datetime.now())[:-7], "Reading the annotation GTF file ...")
transdict = {} #transcript dictionary, key: transcript id, value: exon bounds 
for line in anno:
    record = line.split()
    if line[0]=="#" or record[0]=='track' or record[2]!='exon': continue
    if record[idloc-1] != "transcript_id": #9-12: gene_id "X";transcript_id "Y";
        print("transcript_id not in expected location, exiting"); sys.exit()
    transcript = record[idloc].strip('";')
    coord = [record[0], record[6], int(record[3]), int(record[4])]
    if transcript not in transdict: transdict[transcript] = [coord]
    else: transdict[transcript].append(coord)
#print(transdict)


#############Store all the junctions in a dictionary, each junction as a tuple
#allow shifts to the left or right, e.g. by 2nt (set by the shift option). 
print(str(datetime.now())[:-7], "Building the splicing junction database ...")
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


print(str(datetime.now())[:-7], "Started filtering splice junctions ...")
#############Find locations of all gaps, only N
#tolerate ambiguities at the junctions if off by 2nt as defined by shift

#for line in inputsam:
#    if line[0] == "@": outstring += line; continue

for align in inputbam:
    #QNAME = align.query_name  #line[0]
    CIGAR = align.cigarstring  #line[5]
    #FLAG = align.flag
    #SEQ   = align.query_sequence
    #QUAL  = align.qual
    
    line = align.to_string()
    #if line[0] == "@": outstring += line; continue

    inputcount +=1
    if not inputcount%1000000:
        print(str(datetime.now())[:-7], "Processed", inputcount, "reads ...")
        #outputsam.write(outstring); outstring='' #write output to free up memory
    
    #if 'N' not in line.split()[5]: continue
    if 'N' not in CIGAR: continue
    gaps = getgaps(line); gapcount +=1

    #short gaps are tolerated.
    goodidex = [0]
    if not short: #only remove alignments that only have splice junctions.
        goodidx=[0 if gap in junctdict else 1 for gap in gaps]
    elif short: #remove alignments where all gaps are short or spliced
        goodidx=[0 if gap in junctdict or gap[3]-gap[2]<4 else 1 for gap in gaps]
    if max(goodidx)==1: #0: bad gaps (spliced or 1-2nt); 1: good gaps 
        #at least one gap is good. 
        #outstring+=line;
        outputbam.write(align);
        ligcount+=1 
        if len([i for i in goodidx if i==1])>1: ligmcount+=1
        
anno.close()
#inputsam.close()
#outputsam.write(outstring)
#outputsam.close()
inputbam.close()
outputbam.close()
#logfile.close()

                   
print("\n                        Total alignments:", inputcount)
print("                   All gapped alignments:", gapcount)
print("    Alignments with at least 1 good gaps:", ligcount)
print("    Alignments with at least 2 good gaps:", ligmcount)
print("  Number of annotated splicing junctions:", len(junctdict)/5, '\n')
print(str(datetime.now())[:-7],"Finished gapfilter.py successfully\n")

"""
Analysis of AMT_Stress_trim_nodup_bc07_starhg38Aligned_prim_N.nosplice.sam
                        Total alignments: 3570138
                   All gapped alignments: 3570138
Alignments with at least 1 nonsplice gap: 3549243
  Number of annotated splicing junctions: 261919
It seems that spliced reads are not a major problem in my final output. 
"""
