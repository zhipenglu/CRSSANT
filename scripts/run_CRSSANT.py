#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# RNA crosslinking, proximity ligation and high throughput sequencing produces non-continuous reads that indicate base pairing and higher order interactions,\
# either in RNA secondary structures or intermolecular complexes. CRSSANT (pronounced 'croissant') is a computational pipeline for analyzing non-continuous/gapped \
# reads from a variety of methods that employ the crosslink-ligation principle, including PARIS, LIGR, SPLASH, COMRADES, hiCLIP, etc. CRSSANT optimizes short-read \
# mapping, automates alignment processing, and clusters gap1 and trans alignments into duplex groups (DG) and non-overlapping groups (NG). More complex arrangments \
# are assembled into higher level structures. In particular gapm alignments with 2 gaps or 3 segments are assembled into tri-segment groups (TGs). Overlapping alignments \
# are used to discover homotypic interactions (RNA homodimers).

# Briefly, the CRSSANT pipeline operates as follows. First, sequencing reads that have been processed to remove adapters are mapped references with STAR and a new \
# set of optimized options. Second, alignments are filtered, rearranged and classified into different types (gaptypes.py and gapfilter.py). Third, we use network \
# analysis methods to cluster non-continuous alignments into DGs and calculate the confidence for each DG. The DGs are used as the foundation for the assembly of TGs.

# CRSSANT is written in Python and available as source code that you can download and run directly on your own machine (no compiling needed). An earlier version of \
# the DG assembly method is available here: (https://github.com/ihwang/CRSSANT). For more about the CRSSANT pipeline, please see the bioRxiv preprint by Zhang et al. 2021.

"""
Created on Tue Nov 01 2021

@author: Minjie Zhang (minjiez@usc.edu)

Wrapper script for the 2 parts of the CRSSANT pipeline:
Step 1: Classsify alignments
Step 2: Segment and gap statistics
Step 3: Filter spliced and short gaps
Step 4: Cluster gap1 and trans alignments to DGs
Step 5: Cluster gapm alignments to TGs
If bam files are provided, the software will convert bam to sam fisrt.

This version of the pipeline is designed to be run on a local machine.
"""

import os
import argparse
import time
import sys
PATH=(sys.path[0])+'/scripts'


def Run_CRSSANT():
    start_time = time.time()

    bam2sam_cmd = 'samtools view -h {} > {}'
    bam2sam_cmd = bam2sam_cmd.format(input_dir+sample_name+'.bam', input_dir+sample_name+'.sam')
    
    gaptype_cmd = 'python {}/gap_types.py {} {} {} {} {}'
    gaptype_cmd = gaptype_cmd.format(PATH, input_dir+sample_name+'.sam', output_dir_types+outprefix, glenlog, minseglen, npro)

    gap1filter_cmd = 'python {}/gap_filter.py {} {} {} {} {} gap1'
    gap1filter_cmd = gap1filter_cmd.format(PATH, Gtf, output_dir_types+outprefix+'gap1.sam',
                                         output_dir_types+outprefix+'gap1_filtered.sam', idloc, re_shortGap)
    
    gapmfilter_cmd = 'python {}/gap_filter.py {} {} {} {} {} gapm'
    gapmfilter_cmd = gapmfilter_cmd.format(PATH, Gtf, output_dir_types+outprefix+'gapm.sam',
                                         output_dir_types+outprefix+'gapm_filtered.sam', idloc, re_shortGap)
    
    gaplendist_cmd = 'python {}/gap_lendist.py {} {} all'
    gaplendist_cmd = gaplendist_cmd.format(PATH, output_dir_types+outprefix+'gap1.sam', output_dir_statistics+outprefix)
    
    seglendist_cmd = 'python {}/seg_lendist.py {} {}'
    seglendist_cmd = seglendist_cmd.format(PATH, output_dir_types+outprefix+'gap1.sam', output_dir_statistics+outprefix)
    
    crssant_pre_cmd = 'python {}/crssant_pre.py {} {} {}'
    crssant_pre_cmd = crssant_pre_cmd.format(PATH, output_dir_types+outprefix+'gap1_filtered.sam',
                                             output_dir_types+outprefix+'trans.sam', output_dir_DGs+outprefix)
    
    crssant_cmd = 'python {}/crssant.py -cluster {} -t_o {} -out {} {} {} {}'
    crssant_cmd = crssant_cmd.format(PATH, cluster, t_o, output_dir_DGs, output_dir_DGs+outprefix+'_crssant.sam', genesfile,
                                     output_dir_DGs+outprefix+'_crssant_plus.bedgraph,'+output_dir_DGs+outprefix+'_crssant_minus.bedgraph')
    
    gapm_cluster_cmd = 'python {}/gapm_cluster.py {} {} {}'
    gapm_cluster_cmd = gapm_cluster_cmd.format(PATH, output_dir_DGs+outprefix+'_crssant.'+cluster+'.t_o'+str(t_o)+'_dg.bedpe',
                                               output_dir_types+outprefix+'gapm_filtered.sam', output_dir_TGs+outprefix)
    
    # Check if FASTQ option was specified. If so, run mapping
    if not os.path.exists(input_dir+sample_name+'.sam'):    os.system(bam2sam_cmd)
    os.system(gaptype_cmd)
    os.system(gap1filter_cmd)
    os.system(gapmfilter_cmd)
    os.system(gaplendist_cmd)
    os.system(seglendist_cmd)
    os.system(crssant_pre_cmd)
    os.system(crssant_cmd)
    os.system(gapm_cluster_cmd)

    end_time = time.time()
    time_taken = round((end_time - start_time) / 60, 2)
    print('Total time taken: {} mins'.format(time_taken))


if __name__ == '__main__':
    # Inputs from user
    parser = argparse.ArgumentParser(description='Run CRSSANT')
    parser.add_argument('input_dir', help='Path to input directory')
    parser.add_argument('output_dir', help='Path to output directory')
    parser.add_argument('sample_name', help='Name of sample')
    parser.add_argument('Gtf', help='Name of Gtf annotation file')
    parser.add_argument('idloc', help='Column num of the transcript_id (space sep.) in Gtf file')
    parser.add_argument('genesfile', help='Name of gene bed file')
    parser.add_argument('outprefix', help='outprefix')

    args = parser.parse_args()
    input_dir = args.input_dir
    output_dir = args.output_dir
    sample_name = args.sample_name
    Gtf = args.Gtf
    idloc = int(args.idloc)
    genesfile = args.genesfile
    outprefix = args.outprefix

    # ----- Inputs not specified by the user. Modify these as needed. ------- #
    # Inputs for Step 1 - Classsify alignments
    glenlog = -1 # Scaling factor for gap extension penalty, equivalent to scoreGenomicLengthLog2scale in STAR
    minseglen = 15  # Minimal length for a segment in an alignment to be considered confident for building the connection database
    npro = 8 # Number of CPUs to use for the run, depending availability of resources.
    
    # Inputs for Step 2 - Filter spliced and short gaps
    #idloc = 11 # column num of the transcript_id (space sep.) in Gtf file'
    re_shortGap = 'yes' # 'yes' to remove alignments with only 1/2nt gaps, or 'no'
    
    # Inputs for Step 3 - DG assembly
    cluster = 'cliques' # clustering method, "cliques" or "spectral"
    t_o = 0.1 # Overlap threshold 0-1 inclusive. Default: 0.5 for "spectral", 0.1 for "cliques"
    t_eig = 5 # Eigenratio threshold (positive) for "spectral" only. Default: 5
    
    # Make sure the 'input' directory (inside the 'project' directory)
    # contains the appropriate files inside them.
    # It is not necessary to manually create the 'output' directory.
    if not input_dir.endswith('/'): input_dir += '/'
    if not output_dir.endswith('/'): output_dir += '/'
    
    output_dir_types = output_dir + 'alignments_classsify' + '/'
    output_dir_DGs = output_dir + 'alignments_DGs'+ '/'
    output_dir_TGs = output_dir + 'alignments_TGs'+ '/'
    output_dir_statistics = output_dir + 'alignments_statistics'+ '/'
    if not os.path.exists(output_dir):  os.makedirs(output_dir)
    if not os.path.exists(output_dir_types):  os.makedirs(output_dir_types)
    if not os.path.exists(output_dir_DGs):  os.makedirs(output_dir_DGs)
    if not os.path.exists(output_dir_TGs):  os.makedirs(output_dir_TGs)
    if not os.path.exists(output_dir_statistics):  os.makedirs(output_dir_statistics)

    Run_CRSSANT()
