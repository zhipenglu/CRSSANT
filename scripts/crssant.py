"""
CRSSANT: Crosslinked RNA Secondary Structure Analysis using Network Techniques
Main script for running CRSSANT analysis and discovery pipelines.
Author: Irena Fischer-Hwang and Zhipeng Lu
Contact: ihwang@stanford.edu, zhipengl@usc.edu.
Finalized on 2020. 
Currently only works under python 3.6 because networkx is only
compatible with python 3.6, but not newer versions. The current networkx version
2.4 is not compatible with current python 3.8. 

list of all functions in this file in 6 sections roughly in logical order.
1. Preprocessing: create dictionaries for reads, genes and alignments
2. Subfunctions: process cigar, calculate overlap and get gene information.  
3. Graphing: create a weighted graph, spectral clustering and get cliques 
4. DG analysis: reorganize and filter the duplex groups
5. Output: set up the output format. 
6. Parse arguments and set up multiprocessing run. 

Preprocessing, convert input to dictionaries. 5 functions
def getreads(alignfile):
    Function to create dictionary k=QNAME, v=[alignments]
def getgenes(genesfile): 
    Function to create dictionary k=(chr,strand,pos), v=[names]
def line2info(line):
    Function to convert 1 line from a chimera to a list of information
def getalign(readsdict,genesdict):
    Function to create dictionary k=(gene1,gene2), v=[[alignid,info1,info2]...]
def getcov(alignments):
    Function to calculate coverage across the genome k=(chr,strand), v=[values]

Subfunctions. 3 functions
def get_overlaps(inds_1, inds_2):
    Calculate the overlaps between 2 alignments and the distance spanned by them
def get_overlap_ratios(inds_1, inds_2):
    Function to calculate overlaps ratio between 2 alignments: ratio_l, ratio_r
def process_cigar(cigar_str):
    Function to process CIGAR string in sequencing reads. return ops, lens
    
Graphing. 3 functions.
def subsample(genealign, covlimit):
    Function to sub sample alignments to speed up graph building and clustering
def graph_align(genealignlimit, t):
    Function that creates a weighted graph representation of alignments: graph
def get_spectral(graph, n=10, t=5):
    Function that performs spectral clustering on graph to reads_dg_dict
def get_cliques(graph):
    Function to gets cliques from subgraphs to reads_dg_dict

DG analysis. 4 functions
def get_dgs(reads_dict, reads_dg_dict):
    Function that creates inverse dictionary of reads_dg_dict: dg_reads_dict
def addextra(dg_align_dict, genealign, geneextra, t_o):
    Function to add the extra alignments to established DGs. 
def filter_dgs(dg_reads_dict, reads_dict):
    Function to filter out invalid DGs. return filtered_dict
def create_dg_dict(dg_reads_dict, reads_dict):
    Function that creates the DG dictionary. return dg_stats_dict
def check_nonoverlapping_reads(dg_reads_dict, reads_dict): #deprecated
    return n_nonoverlap

Output. 2 functions
def write_sam(outsam, dg_align_dict):
    Function that writes the ouput SAM file
def write_dg(dg_file, dg_dict, region):
    Function that writes the output bed format for DG information

Run the analysis. 3 functions
def run_analysis(instance): This function runs the main analysis in parallel.
def parse_args(): Function to parse the arguments
def main(): Main function. Run preprocessing, parallel graph analysis and output

"""

import sys
import argparse
import datetime
import re
import numpy as np
import random
import multiprocessing
import networkx as nx
import scipy as sp
from math import floor, ceil
from itertools import chain, product
from sklearn.cluster import KMeans
#import ushuffle, RNA #not used at the moment



###############################################################################
###############################################################################
#1. Process input. 5 functions: getreads,getgenes,line2info,genealign

def getreads(alignfile): #verified
    """
    function to build a dictionary of all reads. 
    input: gap1.sam and rri.sam from gapanalysis.py and gapfilter.py
    these two types of input are automatically detected. 
    returns readsdict: {readid:[alignment lines], ...}
    """
    readsdict={}
    header = []
    f = open(alignfile, 'r')
    readcount, aligncount = 0, 0
    for line in f:
        if line[0]=='@': header.append(line); continue
        QNAME = line.split('\t')[0]
        if QNAME not in readsdict: readsdict[QNAME]=[]; readcount+=1
        readsdict[QNAME].append(line); aligncount+=1
    print("Number of reads:", readcount)
    print("Number of alignments:", aligncount)
    samheader = ''.join(header)
    return readsdict, samheader

def getgenes(genesfile):
    """
    make a sorted list for genes from bed input
    input BED6: chrom chromStart chromEnd name score strand
    genesdict: (chr,strand,pos): [names] #genes may overlap, pos in 10nt intvls
    """
    genesdict = {} #{(chrom,strand,pos): [names] ...}
    f = open(genesfile, 'r')
    genecount,intvlcount = 0,0
    for line in f:
        record= line.rstrip().split('\t'); genecount+=1
        chrom,name,strand = record[0],record[3],record[5]
        chromStart,chromEnd = int(record[1]),int(record[2])
        for i in range(ceil(chromStart/10)*10, chromEnd, 10): #pos: 10nt intvls
            if (chrom,strand,i) not in genesdict: genesdict[(chrom,strand,i)]=[]
            genesdict[(chrom,strand,i)].append(name); intvlcount+=1
    f.close()
    print("Number of genes:", genecount)
    print("genesdict size:", intvlcount)
    return genesdict    

def line2info(line):
    """convert 1 chimeric alignment line to a list of information"""
    align=line.split('\n')
    RNAME,POS,CIGAR=align[2],int(align[3]),align[5]
    Rlen=sum([int(i[:-1]) for i in re.findall('\d+[MD=X]',CIGAR)])
    STRAND='-' if '{0:012b}'.format(int(align[1]))[-5]=='1' else '+'
    info = [RNAME,STRAND,POS,POS+Rlen,line,[]]
    for i in [floor(POS/10)*10,ceil((POS+Rlen)/10)*10]:
        if (RNAME,STRAND,i) in genesdict:
            info[-1].append(genesdict[(RNAME,STRAND,i)])
    return info #[RNAME,STRAND,POS,POS+Rlen,line,[names]]
        
def getalign(readsdict, genesdict):
    """
    function to build a dictionary of noncontinous alignments
    requires itertools.product, line2info()
    export alignments outside of gene models to another sam file
    input readsdict: {readid:[alignments (str)], ...}
    input genesdict: {(chrom,strand,pos): [name ...] ...}
    aligndict format for each item (line1 and line 2 identical for chimera):
    alignid:[chrom1,strand1,start1,end1,line1,genes1,
             chrom2,strand2,start2,end2,line2,genes2]
    """
    aligndict={}
    for QNAME in readsdict:
        lines = readsdict[QNAME] #all alignments for this read
        if lines[0].split()[-1][:2]=='SA': #SA:Z:RNAME,POS,STRAND,CIGAR,MAPQ,NM;
            alignid = lines[0].split()[0]
            info1,info2 = line2info(lines[0]),line2info(lines[1])
            aligndict[alignid] = info1+info2; continue
        linecount=0 
        for line in lines: #now reads from gap1.sam, normal alignment with 1 gap
            linecount+=1
            align=line.split()
            RNAME,POS,CIGAR=align[2],int(align[3]),align[5]
            STRAND='-'if'{0:012b}'.format(int(align[1]))[-5]=='1' else '+'
            glen=int(re.findall('\d+N', CIGAR)[0][:-1]) #gap length
            segs=[i.rstrip('0123456789') for i in CIGAR.split('N')]
            Rlens=[sum([int(i[:-1]) for i in re.findall('\d+[MD=X]',seg)])
                   for seg in segs]
            POS2=POS+glen+Rlens[0]
            info1 = [RNAME,STRAND,POS,POS+Rlens[0],line,[]]
            info2 = [RNAME,STRAND,POS2,POS2+Rlens[1],line,[]]
            for i in [floor(POS/10)*10,ceil((POS+Rlens[0])/10)*10]:
                if (RNAME,STRAND,i) in genesdict:
                    info1[-1].extend(genesdict[(RNAME,STRAND,i)]) #add genename
            for i in [floor(POS2/10)*10,ceil((POS2+Rlens[1])/10)*10]:
                if (RNAME,STRAND,i) in genesdict:
                    info2[-1].extend(genesdict[(RNAME,STRAND,i)]) #add genename
            alignid = line.split()[0] + '_' + str(linecount)
            aligndict[alignid] = info1+info2
            
    #now build a dictionary of genepairs and alignments: genealigndict
    #(gene1,gene2):[[alignid,info1,info2]...]. info:chrom,strand,start,end,line
    genealigndict = {}
    for key in aligndict:
        align = aligndict[key]
        genes1,genes2 = align[5],align[11]
        if not genes1 or not genes2: continue #export the alignments later
        pairs=list(set([tuple(sorted(list(i)))
                        for i in product(set(genes1),set(genes2))]))
        for pair in pairs:
            if pair not in genealigndict: genealigndict[pair]=[]
            genealigndict[pair].append([key]+align[:5]+align[6:11])
    return genealigndict

def getcov(bedgraphplus,bedgraphminus):
    #Function to make a dictionary of genome coverage in 10nt intervals
    #bedgraph format: chrom, start, end, value
    #Returns covdict: (chrom,strand,pos):maxvalue 
    f,g = open(bedgraphplus, 'r'),open(bedgraphminus, 'r'); covdict={}
    for line in f: #plus strand
        chrom,start,end,value = tuple(line.strip().split())
        range10nt = range(round(int(start)/10),round(int(end)/10))
        for i in [10*i for i in range10nt]:
            if (chrom,'+',i) not in covdict: covdict[(chrom,'+',i)]=0
            covdict[(chrom,'+',i)] = max(covdict[(chrom,'+',i)],int(value))
    for line in g: #minus strand
        chrom,start,end,value = tuple(line.strip().split())
        range10nt = range(round(int(start)/10),round(int(end)/10))
        for i in [10*i for i in range10nt]:
            if (chrom,'-',i) not in covdict: covdict[(chrom,'-',i)]=0
            covdict[(chrom,'-',i)] = max(covdict[(chrom,'-',i)],int(value))
    f.close(); g.close(); return covdict

"""
testcase for the hsrRNA data: 
bedtools genomecov -bg -split -strand + -ibam hsrRNA_reads.bam -g hsrRNA.genome\
>hsrRNA_plus.bedgraph; bedtools genomecov -bg -split -strand - -ibam \
hsrRNA_reads.bam -g hsrRNA.genome >hsrRNA_minus.bedgraph
python 
covdict = getcov('hsrRNA_plus.bedgraph', 'hsrRNA_minus.bedgraph')
f = open('test2.bedgraph', 'w')
for i in sorted(list(covdict.items())):
    f.write('\t'.join([i[0][0],str(i[0][2]-5),str(i[0][2]+5),str(i[1])])+'\n')
"""
###############################################################################
###############################################################################










###############################################################################
###############################################################################
#2. Subfunctions. 3 functions: get_overlaps, get_overlaps_ratios, process_cigar

def get_overlaps(inds_1, inds_2): #verified.
    """
    Calculate the overlaps between two reads and the distance spanned by them.
    The reads are each assumed to comprise two arms, a left and a right arm.
    inds_1 and inds_2: [left start, left stop, right start, right stop]
    Returns overlap_l, overlap_r, span_l, span_r: int, int, int, int
    """
    overlap_l = min(inds_1[1], inds_2[1]) - max(inds_1[0], inds_2[0]) + 1
    overlap_r = min(inds_1[3], inds_2[3]) - max(inds_1[2], inds_2[2]) + 1
    span_l = max(inds_1[1], inds_2[1]) - min(inds_1[0], inds_2[0]) + 1
    span_r = max(inds_1[3], inds_2[3]) - min(inds_1[2], inds_2[2]) + 1
    return overlap_l, overlap_r, span_l, span_r

def get_overlap_ratios(inds_1, inds_2): #verified.
    """
    Function to calculate the overlap ratio between two reads
    inds_1 and inds_2: [left start, left stop, right start, right stop]
    Returns ratio_l, ratio_r : float, float
    The integer divisions are not backward compatible with python2. 
    """
    overlap_l, overlap_r, span_l, span_r = get_overlaps(inds_1, inds_2)
    ratio_l = overlap_l / span_l; ratio_r = overlap_r / span_r
    return ratio_l, ratio_r

def process_cigar(cigar_str): #verified. all CIGAR operations considered. 
    """
    Function to process CIGAR string in sequencing reads
    cigar_str: str. CIGAR string
    Returns: ops, lens : list, list
    All cigar operations are considered. Why set the start and end to 'S'?
    requires re. 
    """
    # Parse the cigar string into operations (ops) and operation lengths (lens)
    ops_raw = re.findall('\D+', cigar_str)
    lens_strs = re.findall('\d+', cigar_str)
    lens_raw = [int(i) for i in lens_strs]
    # Merge duplicate consecutive operations
    ops, lens = [], []; lens.append(lens_raw[0]); ops.append(ops_raw[0])
    for i in range(1, len(ops_raw)): 
        if ops_raw[i] == ops_raw[i-1]: lens[-1] += lens_raw[i]
        else: lens.append(lens_raw[i]); ops.append(ops_raw[i])
    # Standardize all cigar strings to start and end with soft-clipped regions
    if ops[0] != 'S': lens = [0] + lens; ops = ['S'] + ops
    if ops[-1] != 'S': lens = lens + [0]; ops = ops + ['S']
    return ops, lens
###############################################################################
###############################################################################









###############################################################################
###############################################################################
#3. Graphing. 4 functions: subsample, graph_align, get_spectral, get_cliques
#sections 3 and 4 process each (gene1,gene2) individually for parallelization

def subsample(genealign, covlimit): #requires random. 
    #genealign: [[alignid,info1,info2]...]. info:chrom,strand,start,end,line
    cov={}; nalign=len(genealign)
    intvls = [i[1:5] for i in genealign] + [i[6:10] for i in genealign]
    for intvl in intvls:
        for i in intvl[2:4]:
            if (*intvl[:2],i) not in cov: cov[(*intvl[:2],i)]=1
            else: cov[(*intvl[:2],i)]+=1
    maxcov=max(list(cov.values()))
    genealignlimit = genealign #first set the limited list to all alignments
    geneextra = [] #store nonselected alignments. 
    if maxcov>covlimit: 
        random.seed(0)
        idx=random.sample(range(nalign),int(nalign*covlimit/maxcov))
        geneextra = [genealign[i] for i in range(nalign) if i not in idx]
        genealignlimit = [genealign[i] for i in idx]
    return genealignlimit, geneextra


def graph_align(genealign, t): #requires networkx. only graph the sub sampled
    """
    Function that creates a weighted graph representation of the alignments
    Alignments are represented by nodes. Two alignments with left and right
    overlap ratios > threshold t are connected by an edge of weight overlap/span
    input genealign is one item from the genealigndict dictionary.
    (gene1,gene2):[[alignid,info1,info2]...]. info:chrom,strand,start,end,line
    alignids: np array. Read IDs
    genealign: list of alignments
    t: float. Overlap threshold, default 0.1 for cliques and 0.5 for spectral
    Returns graph: NetworkX graph
    requires networkx, get_overlap_ratios(), etc.
    """
    graph = nx.Graph()
    alignids = [i[0] for i in genealign]
    alignindsdict = dict([(i[0],i[3:5]+i[8:10]) for i in genealign])
    #alignindsdict format: (alignid: [start1, end1, start2, end2])
    graph.add_nodes_from(alignids)
    aligninds = np.array(list(alignindsdict.values()))
    sorted_ids = []
    for i in range(4): sorted_ids.append([alignid for (ind, alignid) in
                                          sorted(zip(aligninds[:,i],alignids))])
    for (id_1, inds_1) in zip(alignids, aligninds):
        for i in range(4):
            id_list = sorted_ids[i]; j = id_list.index(id_1); loop_flag = 1
            if i%2 == 0: check_loop = (j < len(aligninds) - 1)
            else: check_loop = (j > 0)
            while (loop_flag == 1) and check_loop:
                if i%2 == 0: j += 1; check_loop = (j < len(aligninds) - 1)
                else: j -= 1; check_loop = (j > 0)
                id_2 = id_list[j]; inds_2 = alignindsdict[id_2]
                if i == 0: check_inds = (inds_2[0] <= inds_1[1])
                elif i == 1: check_inds = (inds_2[1] >= inds_1[0])
                elif i == 2: check_inds = (inds_2[2] <= inds_1[3])
                else: check_inds = (inds_2[3] >= inds_1[2])
                if check_inds:
                    ratio_l, ratio_r = get_overlap_ratios(inds_1, inds_2)
                    if (ratio_l > t) and (ratio_r > t): 
                         graph.add_edge(id_1, id_2, weight=(ratio_l + ratio_r))
                else: loop_flag = 0
    return graph

def get_spectral(graph, n=10, t=5): #This function has bugs I have not fixed. 
    """
    Function that performs spectral clustering on the weighted graph.
    Cluster number, k, is determined by finding the first eigengap that is
    some amount t larger than preceding eigengaps.
    graph : NetworkX graph. Weighted graph representation of all alignments
    n: int. Number of eigenvalues (DG splits) to consider threshold
    t: int. Multiplicity of median eigengap threshold
    Returns align_dg_dict, dg_ind: dict, int
    requires KMeans, networkx functions:
    connected_components(),subgraph.nodes(), subgraph.degree(), etc.  
    """
    dg_ind = 0; align_dg_dict = {}
    subgraphs = [graph.subgraph(c) for c in nx.connected_components(graph)]
    for subgraph in subgraphs:
        k = 1
        if len(subgraph) > 1:
            L=nx.laplacian_matrix(
                subgraph,nodelist=sorted(subgraph.nodes())).todense()
            D=np.diag([subgraph.degree[node]
                       for node in sorted(subgraph.nodes())])
            w, v = sp.linalg.eigh(L, D, type=1)  # Since L always symmetric
            eigengaps = np.diff(w[:(n + 1)])
            if len(eigengaps) > 2:
                if (w[1] > 1) and (w[1] >= 10*np.median(eigengaps[1:])): k = 2
                else:
                    # ignore divide by 0 warning if eigengaps median is 0
                    np.seterr(divide='ignore', invalid='ignore')
                    eigenratios = np.copy(eigengaps)
                    eigenratios[1:] = np.array([
                        eigengaps[i] / np.median(eigengaps[:i])
                        for i in range(1, len(eigengaps))])
                    if max(eigenratios) >= t: k = np.argmax(eigenratios>=t)+2
            Y = np.transpose(v[:k])
            kmeans = KMeans(n_clusters=k, random_state=0).fit(Y)
            kmeans.labels_ += dg_ind
            subgraph_dict = dict(zip(sorted(subgraph.nodes()), kmeans.labels_))
        else: subgraph_dict = {list(subgraph)[0]: dg_ind}
        align_dg_dict = {**align_dg_dict, **subgraph_dict}; dg_ind += k
    return align_dg_dict #asignment of alignment IDs to DG numbers.

def get_cliques(graph):
    """
    Function that gets cliques from subgraphs of connected components.
    graph : NetworkX graph. Weighted graph representation of all alignments
    Returns align_dg_dict : dict, int
    requires networkx (<2.4) functions:
    connected_components(), subgraph(), find_cliques, etc. 
    """
    dg_ind = 0; align_dg_dict = {}
    subgraphs = [graph.subgraph(c) for c in nx.connected_components(graph)]
    for subgraph in subgraphs:
        if len(subgraph) > 1:
            cliques_kept = []
            cliques_all = list(nx.find_cliques(subgraph))
            cliques_all.sort(key=len)
            while len(cliques_all) > 0:
                cliques_nodes=set(
                    [node for clique in cliques_kept for node in clique])
                clique_test = cliques_all.pop()
                if not set(list(clique_test)).intersection(cliques_nodes):
                    cliques_kept.append(clique_test)
            dg_inds=[[dg_ind+i]*len(clique)
                     for i,clique in enumerate(cliques_kept)]
            subgraph_dict = dict(
                zip([node for clique in cliques_kept for node in clique],
                    [ind for inds_list in dg_inds for ind in inds_list]))
            align_dg_dict = {**align_dg_dict, **subgraph_dict}
            dg_ind += len(cliques_kept)
        else: read_id=list(subgraph)[0];align_dg_dict[read_id]=dg_ind;dg_ind+=1
    return align_dg_dict #asignment of alignment IDs to DG numbers. 

###############################################################################
###############################################################################












###############################################################################
###############################################################################
#4. DG analysis. 5 functions: get_dgs, addextra, filter_dgs, create_stats and
#check_nonoverlap. these align/dg dictionaries are very confusing

def get_dgs(align_dg_dict):
    """
    Function that creates inverse dictionary of align_dg_dict
    align_dg_dict: dict. Dictionary of alignments and clustering DG assignments
    Returns dg_align_dict: dict, k=dg_id, v=[alignids]
    align_dg_dict comes from get_spectral(graph) or get_cliques(graph)
    """
    dgs_list = set(align_dg_dict.values()) #list of all duplex groups
    dg_align_dict = {}
    for dg in dgs_list:
        dg_align_list =[x for (x,y) in align_dg_dict.items() if y == dg]
        dg_align_dict[dg] = dg_align_list
    return dg_align_dict
    #test case: 

def addextra(dg_align_dict, genealign, geneextra, t_o):
    #dg_align_dict format: k=dg_id; v=[alignids]
    #requires get_overlaps and get_overlap_ratios
    #genealign: [[alignid,info1,info2]...]. info:chrom,strand,start,end,line
    #get the medians of each DG terminals
    indsdict = dict([(i[0],(i[3:5],i[8:10])) for i in genealign])
    medians={} #format: k=dg_id, v=[start1_m,end1_m,start2_m,end2_m]
    for (dg, alignids) in dg_align_dict.items(): #first calculate medians
        inds = np.array([indsdict[i] for i in alignids])
        dg_inds = sum([list(i) for i in np.median(inds, axis=0)],[])
        medians[dg] = dg_inds
    for align in geneextra:
        overlapdict = {} #overlap with dg medians: k=dg_id, v=(ratio_l,ratio_r)
        for (dg,med) in medians.items(): #add geneextra to the duplex groups
            ratio_l,ratio_r = get_overlap_ratios(align[3:5]+align[8:10],med)
            if ratio_l>t_o and ratio_r>t_o: overlapdict[dg] = ratio_l*ratio_r
        if overlapdict:
            dg_align_dict[max(overlapdict,key=overlapdict.get)].append(align[0])
    return dg_align_dict #updated dg_align_dict with all clustered alignments

def filter_dgs(dg_align_dict, genealign):
    """
    Function to filter out invalid DGs
    DGs whose alignments are all identical are eliminated.
    Here we compare (chr1,strand1,start1,end1,chr2,strand2,start2,end2)
    DGs with fewer than two alignments are also eliminated.
    genealign, format as follows
    (gene1,gene2):[[alignid,info1,info2]...]. info:chrom,strand,start,end,line
    dg_align_dict: Dictionary of DGs and their alignments
    Returns filtered_dict : dict
    """
    filtered_dict = {}
    for (dg, alignids) in dg_align_dict.items():
        num_align = len(alignids)
        infolist = [tuple(i[1:5]+i[6:10]) for i in genealign]
        if num_align>1 and len(set(infolist))>1: filtered_dict[dg] = alignids
    return filtered_dict

def create_stats(dg_align_dict, genealign, covdict):
    """
    Function that creates the DG stats dictionary
    The DG stats dictionary no longer has individual alignments but instead
    contains metadata for each DG, including number of alignments, coverage, and
    nonoverlapping group (NG). Coverage fraction is defined as c/sqrt(a*b) where
        + c = number of alignments in a given DG
        + a = number of alignments overlapping the left arm of the DG
        + b = number of alignments overlapping the right arm of the DG
    dg_align_dict: dict. Dictionary of DGs and their alignments
    genealign: list of alignments to this (gene1,gene2) pair
    format: [[alignid,info1,info2]...]. info:chrom,strand,start,end,line
    Returns dg_stats_dict: dict
    """
    ng_ind = 0
    gene_align = list(chain.from_iterable(dg_align_dict.values()))#all alignids
    dg_stats_dict = {}; ng_dict = {}
    for (dg, alignids) in dg_align_dict.items():
        # Get DG indices: dg_min, dg_max and dg_inds (4 medians for each dg)
        dg1aligndict = dict.fromkeys(alignids)
        dg_min=min([i[3] for i in genealign if i[0] in dg1aligndict])
        dg_max=max([i[9] for i in genealign if i[0] in dg1aligndict])
        inds=np.array([i[3:5]+i[8:10]for i in genealign if i[0]in dg1aligndict])
        dg_inds = [int(i) for i in np.median(inds, axis=0)]

        #Calculate coverage based on non-continuous alignments gap1+rri
        chrom1,strand1,chrom2,strand2=tuple(genealign[0][1:3]+genealign[0][6:8])
        start1,end1,start2,end2=tuple(dg_inds)
        range1=[i*10 for i in range(floor(int(start1)/10),ceil(int(end1)/10))]
        range2=[i*10 for i in range(floor(int(start2)/10),ceil(int(end2)/10))]
        overlap1=max([covdict[(chrom1,strand1,i)] for i in range1 if
                      (chrom1,strand1,i) in covdict])
        overlap2=max([covdict[(chrom2,strand2,i)] for i in range2 if
                      (chrom2,strand2,i) in covdict])
        covfrac=len(alignids)/np.sqrt(overlap1*overlap2)
        #print(dg,len(alignids),overlap1,overlap2,covfrac)

       
        # Assign NG for compact visualization in genome browsers e.g. IGV
        # ng_dict format: ng_id: [dg_id list]
        # ng_dg_align format: list of DG ids.
        # ng_dgs: list of dgs for each ng.
        # DGs in each (gene1,gene2) pair are processed separately.
        # To visualize alignments properly in IGV,
        # set "Sort alignments by tag DG" and "Group alignments by tag NG"
        if len(ng_dict)==0: ng_dict[ng_ind]=[dg]; dg_ng=ng_ind; ng_ind+=1
        else:
            ng_assigned = 0
            for (ng, ng_dgs) in ng_dict.items():
                ng_overlaps = np.zeros(len(ng_dgs))
                for (i, ng_dg) in enumerate(ng_dgs):
                    ng_dg_align = dict.fromkeys(dg_align_dict[ng_dg])
                    dginfo=[i for i in genealign if i[0] in ng_dg_align]
                    #dginfo: [[alignid,info1,info2]...].
                    #info:chrom,strand,start,end,line
                    ng_dg_min = min([i[3] for i in dginfo])
                    ng_dg_max = max([i[9] for i in dginfo])
                    overlap = min(dg_max, ng_dg_max) - max(dg_min, ng_dg_min)
                    if overlap > 0: ng_overlaps[i] = 1
                if sum(ng_overlaps) == 0:
                    ng_dict[ng].append(dg); dg_ng = ng; ng_assigned += 1; break
            if ng_assigned == 0:
                ng_dict[ng_ind] = [dg]; dg_ng = ng_ind; ng_ind += 1
        
        dg_stats_dict[dg]={'chromstrand': [chrom1,strand1,chrom2,strand2],
                           'dg_inds': dg_inds, 'covfrac': covfrac,
                           'num_align': len(alignids), 'NG': dg_ng}
    return dg_stats_dict

###############################################################################
###############################################################################











###############################################################################
###############################################################################
#5. Output. 2 functions: write_sam, write_dg.

def writesam(outsam,genealigndict,dg_filtered_dict,dg_stats_dict,genepair):
    """
    Function that writes the ouput SAM file
    SAM file includes only reads that contributed to DG groups, i.e. if DG had
    only one read, which does not pass the DG filter and the read is not output
    genealigndict format for each item: 
    #(gene1,gene2):[[alignid,info1,info2]...]. info:chrom,strand,start,end,line
    dg_filtered_dict: Dictionary of DG and alignments {dg_id:align_id}
    dg_stats_dict: Dictionary of DG metadata, arm_inds,covfraction,num_align,NG
    genepair: tuple, (gene1,gene2). Intramolecular structure when gene1==gene2
    pseudocode: for each dg in dg_filtered_dict, get alignments from aligndict,
    and ng information from dg_stats_dict, and then write alignments to outsam
    """
    readsdict = {}
    f= open(outsam, 'a')
    aligndict = dict([(i[0],(i[5],i[10])) for i in genealigndict[genepair]])
    for (dg,alignids) in dg_filtered_dict.items():
        ng = dg_stats_dict[dg]['NG']
        dgngtags = '\tDG:Z:'+','.join(genepair)+','+str(dg)+'\tNG:i:'+str(ng)
        align = [] #all alignments for this DG
        for alignid in alignids:
            align.append(aligndict[alignid][0].strip()+dgngtags)
            if aligndict[alignid][0] !=aligndict[alignid][1]:
                align.append(aligndict[alignid][1].strip()+dgngtags)
        f.write('\n'.join(align)+'\n')
    f.close()
    return


def writedg(outdg,dg_stats_dict,genepair): #bedpe format output
    """
    see bedpe definition here: 
    https://bedtools.readthedocs.io/en/latest/content/general-usage.html
    chrom1 start1 end1 chrom2 start2 end2 name score strand1 strand2 additional
    bedpe can be easily converted to the bed12 format to visualize arcs. 
    """
    with open(outdg, 'a') as f:
        for (dg, dg_info) in dg_stats_dict.items():
            chromstrand,dg_inds = dg_info['chromstrand'],dg_info['dg_inds']
            covfrac,num_align = dg_info['covfrac'], dg_info['num_align']
            dg_str = '{0},{1},{2:.3f}'.format(','.join(genepair), dg, covfrac)
            line = [chromstrand[0], str(dg_inds[0]), str(dg_inds[1]),
                    chromstrand[2], str(dg_inds[2]), str(dg_inds[3]), dg_str, 
                    str(num_align), chromstrand[1], chromstrand[3]]
            f.write('\t'.join(line)+'\n')
    return
###############################################################################
###############################################################################










###############################################################################
###############################################################################
#6. Run the analysis. 3 Functions: parse_args, run_analysis and main

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser=argparse.ArgumentParser(
        description='CRSSANT groups PARIS reads into DGs and SGs',
        epilog='python ./CRSSANT.py align.sam ref.bed ref.fa out')
    parser.add_argument('alignfile', help='Path to alignedment file (SAM)')
    parser.add_argument('genesfile', help='Path to gene annotations (BED)')
    parser.add_argument('bedgraphs', help='Path to genome coverge (bedgraph)'
                        'Coma separated list of files for the + and - strands')
    parser.add_argument('-out', help='Path of output')
    parser.add_argument('-cluster', help='"cliques" or "spectral"(default)')
    parser.add_argument('-n', help='Number of threads. Default is 8')
    parser.add_argument('-covlimit', help='Max coverage to be directly graphed'
                        'Default: 1000. Extra alignments are added afterwards')
    parser.add_argument('-t_o', help='Overlap threshold 0-1 inclusive'
                        'Default: 0.5 for "spectral", and 0.1 for "cliques"')
    parser.add_argument('-t_eig', help='Eigenratio threshold (positive number)'
                        'Default: 5 for "spectral". Not needed for "cliques"')
    args = parser.parse_args(sys.argv[1:])
    return args

def run_analysis(instance):#run analysis on 1 group of alignments (gene1,gene2)
    genestart = datetime.datetime.now()
    dg_align_dict = None
    dg_stats_dict = None
    genepair,genealign,covdict,covlimit,cluster,t_eig,t_o = instance
    # Run analysis, build graph, cluster and generate the dg output
    genealignlimit,geneextra = subsample(genealign,covlimit)
    graph = graph_align(genealignlimit, t=t_o)
    if cluster == 'cliques': align_dg_dict = get_cliques(graph)
    else: align_dg_dict = get_spectral(graph, t=t_eig)    
    dg_align_dict = get_dgs(align_dg_dict)
    if geneextra:dg_align_dict=addextra(dg_align_dict,genealign,geneextra,t_o)
    dg_filtered_dict = filter_dgs(dg_align_dict,genealign)
    dg_stats_dict = create_stats(dg_filtered_dict,genealign,covdict) #not yet
    return genepair, dg_filtered_dict, dg_stats_dict
    genestop = datetime.datetime.now()
    genetime = genestop - genestart
    #print(genepair, genetime, len(genealign), len(dg_filtered_dict))
    
def main(): #get it to work on a simple case first, before global optimization
    args = parse_args()
    # 1) Read in bed references and alignments
    genesdict = getgenes(args.genesfile)
    readsdict,samheader = getreads(args.alignfile)
    genealigndict = getalign(readsdict, genesdict)
    covdict = getcov(*(args.bedgraphs.split(',')))
    # 2) Check clustering arguments
    if args.covlimit: args.covlimit = int(args.covlimit)
    else: args.covlimit = 1000
    if args.cluster == 'cliques':
        args.t_o = float(args.t_o) if args.t_o else 0.1
        args.t_eig = None
    elif (not args.cluster) or args.cluster == 'spectral':
        args.cluster = 'spectral'
        args.t_o = float(args.t_o) if args.t_o else 0.5
        args.t_eig = float(args.t_eig) if args.t_eig else 5
    else: sys.exit('ABORTED: Valid clustering methods: "cliques" or "spectral"')

    # 3) Check number of threads
    args.n = int(args.n) if args.n else 8
  
    # 4) Run pipeline. Run each gene pair separately. 
    instances=[(genepair,genealigndict[genepair],covdict,args.covlimit,
                args.cluster,args.t_eig,args.t_o)for genepair in genealigndict]
    pool = multiprocessing.Pool(processes=args.n)
    results = pool.imap_unordered(run_analysis, instances)
    
    # 5) Initalize DG files and output SAM and DG files
    file_base = args.alignfile.split('/')[-1].split('.sam')[0]
    if args.cluster=='cliques':
        clustering_str='%s.t_o%s'%(args.cluster,args.t_o)
    else: clustering_str='%s.t_o%s.t_eig%s'%(args.cluster,args.t_o,args.t_eig)
    file_base += '.%s' %clustering_str
    outsam = args.out + file_base + '.sam'
    outdg = args.out + file_base +  '_dg.bedpe'
    with open(outsam, 'w') as f: f.write(samheader)
    with open(outdg, 'w') as f: pass
    for (genepair,dg_filtered_dict,dg_stats_dict) in results:
        writesam(outsam,genealigndict,dg_filtered_dict,dg_stats_dict,genepair)
        writedg(outdg,dg_stats_dict,genepair)
        #write sg output later. DG and NG names are not unique

if __name__ == '__main__':
    multiprocessing.freeze_support()
    main()
sys.exit()
###############################################################################
###############################################################################




