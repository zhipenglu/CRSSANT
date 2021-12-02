"""
minjiezhang123@gmail.com    2020-12-16      python 3

"""


#1. input and output setup
################################################################################
################################################################################
#this section sets up the input and output
import sys, argparse, numpy, os, re, itertools, random, math
from datetime import datetime
import seaborn as sns
import matplotlib.pyplot as plt
from math import floor, ceil
from itertools import product
import matplotlib as mpl

if len(sys.argv) < 3:
    print("Usage:           python homosam_heatmap.py  homo_sam  target_pos  genome.fai  winbin  interval_maker  overlap_cutoff  outprefix")
    print("homo_sam:        homo sam (overlapped chimeras) file generated by gaptypes.py")
    print("gene_bed_file:   gene bed filr")
    print("                 eg: chr1  11  35  ACTB  1000  +")
    print("target_pos:      target position, eg: chr1:1233-1347:+")
    print("genome.fai:      genome length file")
    print("winbin:          winbin for heatmap")
    print("interval_maker:  interval maker (nt) to plot the vlines and hlines")
    print("overlap_cutoff:  cutoff value for overlap length, only the homodimers with overlap regions more than cutoff value will be analyzed")
    print("outputprefix                       ")
    sys.exit()

inputsam = open(sys.argv[1],'r')
target_pos = sys.argv[2]
genome_fai = sys.argv[3]
winbin = int(sys.argv[4])
interval = int(sys.argv[5])
overlap_cutoff=int(sys.argv[6])
outprefix = sys.argv[7]
################################################################################


#2. subfunctions
################################################################################
def timenow(): return str(datetime.now())[:-7]

def mergeCIGAR(CIGAR): 
    #merge all operations that consume the reference, i.e. MI=X
    #example: 1S2M3N4M5I6M7S -> 1S2M3N10M7S 
    ops = re.findall('\d+[MNISH=X]', CIGAR) #all that consume query
    newops = [ops[0]]
    for op in ops[1:]: #concatenate all internal ops that consume query [MIS=X
        if op[-1] not in "I=X":
            if newops[-1][-1]=="M" and op[-1]=="M":
                newops[-1] = str(int(newops[-1][:-1])+int(op[:-1]))+"M"
            else: newops.append(op)
    newCIGAR = ("".join(str(i) for i in newops))
    return newCIGAR

# get homo chimeras segs
def get_homo_segs(line): #turn an alignment into segs [(RNAME, STRAND L, R), ...]
    #A00208:104:HJCGWDRXX:1:2128:6696:28557-HEK293toRNA      2064    chr1    1109945 255     24M2I1D2I1D33M  *       0       0\
    #GGCCGTGAAATGGGGGGTCAGTTAGTGTGCCTGGCAGCGAGCCATCACTGTTTCCATACCG   :FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFF:FFF   NH:i:1  HI:i:1  AS:i:25 nM:i:0  NM:i:0  MD:Z:26 jM:B:c,-1jI:B:i,-1       ch:A:1  SA:Z:chr1,1109969,-,21S35M27S,255,0;
    align = line.rstrip('\n').split('\t')
    readID, CHR, POS, CIGAR = align[0], align[2], int(align[3]), align[5]
    STRAND = '-' if '{0:012b}'.format(int(align[1]))[-5] == '1' else '+'
    seglens = [int(i[:-1]) for i in re.findall('[0-9]+M', CIGAR)] #seg lengths
    overlap = len(re.findall('2I1D', CIGAR))
    seg1 = [CHR, POS, POS+seglens[0]+overlap-1, STRAND, seglens[0]+overlap]
    seg2 = [CHR, POS+seglens[0], POS+seglens[0]+seglens[1]+overlap-1, STRAND, seglens[1]+overlap]
    overlapseg = [CHR, POS+seglens[0], POS+seglens[0]+overlap-1, STRAND, overlap]
    segs = [seg1, seg2, overlapseg]
    return readID,segs

# get gap1 segs
def get_gap1_segs(line): #turn an alignment into segs [(RNAME, STRAND L, R), ...]
    align = line.split(); segs = []; mingap_len = 0
    CHR, POS, CIGAR = align[2], int(align[3]), align[5]
    STRAND = '-' if '{0:012b}'.format(int(align[1]))[-5] == '1' else '+'
    seglens = [int(i[:-1]) for i in re.findall('[0-9]+M', mergeCIGAR(CIGAR))] #seg lengths
    gaplens = [int(i[:-1]) for i in re.findall('[0-9]+N', mergeCIGAR(CIGAR))] #gap lengths
    Ns =[i.rstrip('0123456789') for i in mergeCIGAR(CIGAR).split('M')]
    rx = [] #reference consumed: MD=X
    for N in Ns:
        rx.append(sum([int(i[:-1]) for i in re.findall('[0-9]+[ND=X]', N)]))
    for i in range(len(seglens)): #combine ref and segment lengths to make the junctions
        l, r = POS+sum(rx[:i+1])+sum(seglens[:i]), POS+sum(rx[:i+1])+sum(seglens[:i+1])-1
        if (CHR, int(l), int(r), STRAND) not in segs: segs.append((CHR, int(l), int(r), STRAND))
        segs.sort()
    #print(segs[align[0]])
    mingap_len = min(gaplens)
    return mingap_len,segs

def DG_shuffle(segsdict, genome_fai, chr, start, end, maxcount):
    # output newDGspan in to bed1 and bed2
    segs_shuffled = {}; count = 0
    bed1 = open('bed1.bed','w')
    bed2 = open('bed2.bed','w')
    for readID in segsdict:
        seg_l, seg_r, seg_over = segsdict[readID]
        left_arm = str(seg_l[0])+'\t'+str(seg_l[1])+'\t'+str(seg_l[2])+'\n'
        right_arm = str(seg_r[0])+'\t'+str(seg_r[1])+'\t'+str(seg_r[2])+'\t'+str(readID)+'\n'
        bed1.write(left_arm)
        bed2.write(right_arm)
    bed1.close()
    bed2.close()
    # output incl.bed and genome.bed
    inclbed = open('incl.bed','w')
    inclbed.write(chr+'\t'+str(start)+'\t'+str(end)+'\n')
    inclbed.close()
    # bedtools shuffle
    os.system("bedtools shuffle -incl incl.bed -i bed1.bed -g %s > bed1_shuffled.bed" %(genome_fai))
    os.system("bedtools shuffle -incl incl.bed -i bed1.bed -g %s >> bed1_shuffled.bed" %(genome_fai))
    os.system("bedtools shuffle -incl incl.bed -i bed2.bed -g %s > bed2_shuffled.bed" %(genome_fai))
    os.system("bedtools shuffle -incl incl.bed -i bed2.bed -g %s >> bed2_shuffled.bed" %(genome_fai))
    os.system("paste bed1_shuffled.bed  bed2_shuffled.bed > bed1_2_shuffled.bed")
    # import shuffled
    inputbedpe=open('bed1_2_shuffled.bed','r')
    for line in inputbedpe:
        count += 1
        align = line.rstrip('\n').split('\t')
        keyID = align[6]+'-'+str(count)
        if count<=maxcount:
            segs_shuffled[keyID] = [[align[0],int(align[1]),int(align[2])],[align[3],int(align[4]),int(align[5])]]
    inputbedpe.close()
    os.system("rm incl.bed  bed1.bed  bed2.bed  bed1_shuffled.bed  bed2_shuffled.bed  bed1_2_shuffled.bed")
    return segs_shuffled

# a: [chr1, 1233, 1347, +]
# b: [chr1, 1230, 1377, +]
def getOverlap(a, b):
    seganno = "False"
    if a[0]==b[0] and a[3]==b[3]:
        overlap = int(min(int(a[2]),int(b[2]))) - int(max(int(a[1]),int(b[1])))
        if overlap >1 and overlap >= readslen: seganno = "Ture"
    return seganno

## test whether bed1 was included by bed2:
def bedinregions(bed, target_region):
    # bed = ['hs45S', '3998', '4012', '+']
    # regions = [['hs45S', '3654', '5523', '+'], ['hs45S', '6600', '6757', '+'], ['hs45S', '7924', '12994', '+']]
    bedanno = 'False'
    if str(bed[0]) == str(target_region[0]) and str(bed[3]) == str(target_region[3]):
        if int(bed[1]) >= int(target_region[1]) and int(bed[2]) <= int(target_region[2]):
            bedanno = 'True'
    return bedanno

# get noBP within winbin:
# get_interval(6,12,4,13,'+',5)
# output: [0, 1]
def get_interval(seg_start, seg_end, gene_start, gene_end, gene_strand, winbin):
    start = 1; end = 0; interval = []
    if int(gene_start) <= int(seg_start): start = int(seg_start)-int(gene_start)+1
    if int(seg_end) <= int(gene_end): end = int(seg_end)-int(gene_start)+1
    else: end = int(gene_end)-int(gene_start)+1
    for i in range(floor(start/winbin), ceil(end/winbin), 1):
        if gene_strand == '+':  interval.append(i)
        if gene_strand == '-':  interval.append(ceil((int(gene_end)-int(gene_start)+1)/winbin-1)-i)
    return sorted(interval)

# Python program to illustrate the intersection
# of two lists using set() method
def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))

## Plot heatmap:
cdict = {'red':   [(0.0,  1.0, 1.0),(0.5,  0.8, 0.8),(0.8,  0.5, 0.5),(1.0,  0.0, 0.0)],
         'green': [(0.0,  1.0, 1.0),(0.5,  0.8, 0.8),(0.8,  0.5, 0.5),(1.0,  0.0, 0.0)],
         'blue':  [(0.0,  1.0, 1.0),(0.5,  0.8, 0.8),(0.8,  0.5, 0.5),(1.0,  0.0, 0.0)]}
cmp1 = mpl.colors.LinearSegmentedColormap('name',cdict)

## Plot heatmap:
cdict = {'red':   [(0.0,  1.0, 1.0),(0.4,  0.6, 0.6),(0.6,  0.4, 0.4),(1.0,  0.0, 0.0)],
         'green': [(0.0,  1.0, 1.0),(0.4,  0.6, 0.6),(0.6,  0.4, 0.4),(1.0,  0.0, 0.0)],
         'blue':  [(0.0,  1.0, 1.0),(0.4,  0.6, 0.6),(0.6,  0.4, 0.4),(1.0,  0.0, 0.0)]}
cmp1 = mpl.colors.LinearSegmentedColormap('name',cdict)

def plot_heatmap(heatmap_matrix, heatmap_matrix_shuffle, start, end, interval, winbin, outprefix):
    axes = plt.subplot(2,2,1)
    max_chimeras = int(heatmap_matrix.max())
    #max_chimeras = 4
    sns.heatmap(heatmap_matrix,cmap=cmp1, vmin=0, vmax=max_chimeras, zorder=10)
    for x in range(0, ceil(((int(end)-int(start)+1)/winbin)+(interval/winbin)), int(interval/winbin)): plt.vlines(x, 0, ceil(int(end)-int(start)+1)/winbin, color = "silver", linestyles = "solid", linewidth=1, zorder=10)
    for x in range(0, ceil(((int(end)-int(start)+1)/winbin)+(interval/winbin)), int(interval/winbin)): plt.hlines(x, 0, ceil(int(end)-int(start)+1)/winbin, color = "silver", linestyles = "solid", linewidth=1, zorder=10)
    #plt.plot([0, ceil((int(end)-int(start)+1)/winbin)],[0, ceil((int(end)-int(start)+1)/winbin)], linewidth=2, color = 'silver', zorder=15)
    sns.despine(top=False, right=False, left=False, bottom=False)
    axes.xaxis.set_ticklabels([])
    axes.yaxis.set_ticklabels([])
    plt.xticks([])
    plt.yticks([])
    plt.title(label=outprefix)
    axes = plt.subplot(2,2,2)
    #max_chimeras = ceil(heatmap_matrix_shuffle.max())
    sns.heatmap(heatmap_matrix_shuffle,cmap=cmp1, vmin=0, vmax=max_chimeras, zorder=10)
    for x in range(0, ceil(((int(end)-int(start)+1)/winbin)+(interval/winbin)), int(interval/winbin)): plt.vlines(x, 0, ceil(int(end)-int(start)+1)/winbin, color = "silver", linestyles = "solid", linewidth=1, zorder=10)
    for x in range(0, ceil(((int(end)-int(start)+1)/winbin)+(interval/winbin)), int(interval/winbin)): plt.hlines(x, 0, ceil(int(end)-int(start)+1)/winbin, color = "silver", linestyles = "solid", linewidth=1, zorder=10)
    #plt.plot([0, ceil((int(end)-int(start)+1)/winbin)],[0, ceil((int(end)-int(start)+1)/winbin)], linewidth=2, color = 'silver', zorder=15)
    sns.despine(top=False, right=False, left=False, bottom=False)
    axes.xaxis.set_ticklabels([])
    axes.yaxis.set_ticklabels([])
    plt.xticks([])
    plt.yticks([])
    plt.title(label="shuffled_"+outprefix)
    plt.savefig(outprefix+".pdf")
    plt.show()
    return
################################################################################


#3 start processing the homo sam file.
################################################################################
print(timenow()+" Process homo chimeric sam file...")
target_region = [target_pos.split(':')[0],target_pos.split(':')[1].split('-')[0],target_pos.split(':')[1].split('-')[1],target_pos.split(':')[-1]]
target_reads = {} #{readID: [seg1, seg2, overlapseg]}
overlap_dict_all = {};  overlap_dict_targetregion = {}  #{overlap_len: reads number}


for line in inputsam:
    if line[0]=='@': continue;
    readID,segs = get_homo_segs(line)
    # get the target regions info
    if bedinregions([segs[0][0],segs[0][1],segs[1][2],segs[0][3]], target_region) == 'True':
        if segs[2][4] >= overlap_cutoff:
            # store overlapped region distribution info
            if segs[2][4] not in overlap_dict_targetregion: overlap_dict_targetregion[segs[2][4]] = 0
            overlap_dict_targetregion[segs[2][4]] += 1
            target_reads[readID] = segs
inputsam.close()
#output_bed_cutoff.close()

"""
## convert overlapped bed file to bedgraph file
os.system("sort -k1,1 -k2,2n -o %s %s" % (outprefix+'_morethan2nt.bed', outprefix+'_morethan2nt.bed'))
os.system("bedtools genomecov -i %s -g %s -bg > %s" % (outprefix+'_morethan2nt.bed', genome_fai, outprefix+'_morethan2nt.bedgraph'))
################################################################################


#4 output overlapped length distribution
################################################################################
print(timenow()+" Output overlapped distribution and coverage file...")
output_overlap_dis_targetregion = open(outprefix+'_overlap_distribution_targetregion.txt','w')
for i in sorted(overlap_dict_targetregion):
    output_overlap_dis_targetregion.write(str(i)+'\t'+str(overlap_dict_targetregion[i])+'\n')
output_overlap_dis_targetregion.close()
"""

#5 plot heatmap plots
################################################################################
print(timenow()+" Plot heatmap...")
target_len = int(target_region[2]) - int(target_region[1]) + 1
col_num = ceil(target_len/winbin); row_num = ceil(target_len/winbin)
heatmap_matrix = numpy.zeros((row_num,col_num));            heatmap_matrix_log = numpy.zeros((row_num,col_num))
heatmap_matrix_shuffle = numpy.zeros((row_num,col_num));    heatmap_matrix_log_shuffle = numpy.zeros((row_num,col_num))

# get the heatmap_matrix coverage
for readID in target_reads:
    #print(target_reads[readID])
    seg_l, seg_r, seg_over = target_reads[readID]
    interval_l = get_interval(seg_l[1], seg_l[2], target_region[1], target_region[2], target_region[3], winbin)
    interval_r = get_interval(seg_r[1], seg_r[2], target_region[1], target_region[2], target_region[3], winbin)
    for i in interval_r:
        for j in interval_l:
            heatmap_matrix[i,j] += 1;
            #X: left arm; Y: right arm
    for i in intersection(interval_l, interval_r):
        for j in intersection(interval_l, interval_r):
            heatmap_matrix[i,j] += 1

# get the shuffled heatmap_matrix coverage
target_reads_shuffle = DG_shuffle(target_reads, genome_fai, target_region[0], target_region[1], target_region[2], len(target_reads))
for readID in target_reads_shuffle:
    #print(target_reads[readID])
    seg_l, seg_r = target_reads_shuffle[readID]
    interval_l = get_interval(seg_l[1], seg_l[2], target_region[1], target_region[2], target_region[3], winbin)
    interval_r = get_interval(seg_r[1], seg_r[2], target_region[1], target_region[2], target_region[3], winbin)
    for i in interval_r:
        for j in interval_l:
            heatmap_matrix_shuffle[i,j] += 1;
            #X: left arm; Y: right arm
    for i in intersection(interval_l, interval_r):
        for j in intersection(interval_l, interval_r):
            heatmap_matrix_shuffle[i,j] += 1
                   
# normalize matrix
for i in range(0, row_num, 1):
    for j in range(0, col_num, 1):
        if heatmap_matrix[i,j]==0:  heatmap_matrix_log[i,j]=0
        else: heatmap_matrix_log[i,j] = math.log(heatmap_matrix[i,j],2)
        if heatmap_matrix_shuffle[i,j]==0:  heatmap_matrix_log_shuffle[i,j]=0
        else: heatmap_matrix_log_shuffle[i,j] = math.log(heatmap_matrix_shuffle[i,j],2)
        
# plot heatmap
plot_heatmap(heatmap_matrix_log, heatmap_matrix_log_shuffle, target_region[1], target_region[2], interval, winbin, outprefix)



