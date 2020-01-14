# CRSSANT: Cross-linked RNA Secondary Structure Analysis using Network Techniques
RNA crosslinking, proximity ligation and high throughput sequencing produces non-continuous reads that indicate base pairing and higher order interactions, either in RNA secondary structures or intermolecular complexes. CRSSANT (pronounced 'croissant') is a computational pipeline for analyzing non-continuous/gapped reads from a variety of methods that employ the crosslink-ligation principle, including [PARIS](https://www.ncbi.nlm.nih.gov/pubmed/27180905), [LIGR](https://www.ncbi.nlm.nih.gov/pubmed/27184080), [SPLASH](https://www.ncbi.nlm.nih.gov/pubmed/27184079), [COMRADES](https://www.ncbi.nlm.nih.gov/pubmed/30202058), [hiCLIP](https://www.ncbi.nlm.nih.gov/pubmed/25799984), etc. CRSSANT optimizes short-read mapping, automates alignment processing, and clusters alignments into duplex groups (DG) and non-overlapping groups (NG).

Briefly, the CRSSANT pipeline operates as follows. First, sequencing reads that have been processed to remove adapters are mapped references with STAR and a new set of optimized options. Second, alignments are filtered, rearranged and classified into different types (gaptypes.py and gapfilter.py). Third, we use network analysis methods to cluster non-continuous alignments into DGs and calculate the confidence for each DG.

CRSSANT is written in Python and available as source code that you can download and run on yuor own machine. An earlier version of the DG assembly method is available here: (https://github.com/ihwang/CRSSANT). For more about the CRSSANT pipeline, please see the [bioRxiv preprint by Fischer-Hwang et al.](LINKLINKLINK).



## Table of contents
* [Download and install](https://github.com/zhipenglu/CRSSANT#download-and-install)
* [Step 1: Map reads to the genome](https://github.com/zhipenglu/CRSSANT#map-reads-to-the-genome)
* [Step 2: Classify alignments](https://github.com/zhipenglu/CRSSANT#classify-alignments)
* [Step 3: Filter spliced and short gaps](https://github.com/zhipenglu/CRSSANT#filter-spliced-and-short-gaps)
* [Step 4: Cluster alignments to groups](https://github.com/zhipenglu/CRSSANT#cluster-alignments-to-groups)
* [Test](https://github.com/zhipenglu/CRSSANT#test)


## Download and install
Navigate to the latest [release](https://github.com/zhipenglu/CRSSANT/releases), right click on the source code, and save it to a known path/location, e.g. `CRSSANT_path`. You will need Python version 3.6+ and the following Python packages. We recommend downloading the latest versions of these packages using the Ananconda/Bioconda package manager (follow instructions in links in parentheses):

No special installation is needed, but the dependencies need to be resolved before use. Currently, the NetworkX. Python dependencies are as follows: 
Numpy, SciPy, etc., which are included in Anaconda python 3.6+ NetworkX has not been updated to be compatible with higher python versions. 
* [STAR v2.7.1+](https://github.com/alexdobin/STAR)

* [NetworkX v2.1+](https://networkx.github.io/) ([Anaconda Cloud link](https://anaconda.org/anaconda/networkx))

Additional tools for used for general processing of high throughput sequencing data, including samtools, bedtools, 

* [NumPy](http://www.numpy.org/) ([Anaconda Cloud link](https://anaconda.org/anaconda/numpy))
* [scikit-learn](http://scikit-learn.org/stable/) ([Anaconda Cloud link](https://anaconda.org/anaconda/scikit-learn))
* [SciPy](https://www.scipy.org/) ([Anaconda Cloud link](https://anaconda.org/anaconda/scipy))












## Step 1: Map reads to the genome
It is assumed that the reads have been demultiplexed and adapters removed. Before mapping the reads, genome indices should be generated with the same STAR version. Reads in the fastq format are mapped to the genome using STAR and a set of optimized parameters as follows. `runThreadN` and `genomeLoad` should be adjusted based on available resources and running environment.  

```
STAR --runMode alignReads --genomeDir /path/to/index --readFilesIn /path/to/reads/files --outFileNamePrefix /path/to/output/prefix --runThreadN 1 --genomeLoad NoSharedMemory --outReadsUnmapped Fastx  --outFilterMultimapNmax 10 --outFilterScoreMinOverLread 0 --outSAMattributes All --outSAMtype BAM Unsorted SortedByCoordinate --alignIntronMin 1 --scoreGap 0 --scoreGapNoncan 0 --scoreGapGCAG 0 --scoreGapATAC 0 --scoreGenomicLengthLog2scale -1 --chimOutType WithinBAM HardClip --chimSegmentMin 5 --chimJunctionOverhangMin 5 --chimScoreJunctionNonGTAG 0 -- chimScoreDropMax 80 --chimNonchimScoreDropMin 20
```

Successful STAR mapping generates the following 7 files: `Aligned.out.bam`, `Aligned.sortedByCoord.out.bam`, `Log.final.out`, `Log.out`, `Log.progress.out`, `SJ.out.tab`, and `Unmapped.out.mate1`. The `bam` file is converted back to `sam` for the next step of processing, keeping the header lines (`samtools view -h`). Sorting is not necessary for the next alignment classification step.  

## Step 2: Classify alignments
In this step, alignments in the sam file are filtered to remove low-confidence segments, rearranged and classified into 5 distinct types using `gaptypes.py`. 
```
python gaptypes.py input.sam output_prefix glenlog nprocs
```

Recommended parameters are as follows. 
* `glenlog`: -1. Scaling factor for gap extension penalty, equivalent to `scoreGenomicLengthLog2scale` in STAR 
* `minlen`: 15. Minimal length for a segment in an alignment to be considered confident for building the connection database
* `nprocs`: 10. Number of CPUs to use for the run, depending availability of resources. 

Successful completion of this step results in 7 files. All of these sam files can be converted to sorted bam for visualization on IGV. 
* `cont.sam`: continuous alignments
* `gap1.sam`: non-continuous alignments, each has 1 gap
* `gapm.sam`: non-continuous alignments, each has more than 1 gaps
* `trans.sam`: non-continuous alignments with the 2 arms on different strands or chromosomes
* `homo.sam`: non-continuous alignments with the 2 arms overlapping each other
* `bad.sam`: non-continuous alignments with complex combinations of indels and gaps 
* `log.out`: log file for the run, including input and output alignment counts (for `gap1`, `gapm`, `trans`, `homo` and `bad`)


## Step 3: Filter spliced and short gaps
Output files `gap1.sam` and `gapm.sam` may contain alignments that have only splicing junctions and short 1-2 nt gaps due to artifacts. These are filtered out using `gapfilter.py` before further processing. Splicing junctions and short gaps in other output files can be safely ignored. The `annotation` file containing the splicing junctions should be in GTF format. `idloc`, location of the transcript_ID field, is usually field 11. `short` is set to either `yes` which means 'remove short 1-2nt gaps', or `no`, which means 'ignore short 1-2nt gaps'.  
```
Usage: python gapfilter.py annotation insam outsam idloc short
```
The output from `gap1.sam`, typically named `gap1filter.sam`, only contains alignments that pass the filter. The output from `gapm.sam`, typically named gapmfilter.sam, contain alignments with either 1 or more gaps that pass the filter. The following output is printed to the screen: 

* `Total alignments`
* `All gapped alignments`
* `Alignments with at least 1 good gaps`
* `Alignments with at least 2 good gaps`
* `Number of annotated splicing junctions`

Alignments counts in log.out from step 2 and filtered counts are used to calculate percentage of non-continuous alignments alignments using the following formula: 
```
(gap1filtercount + gapmfiltercount + transcount + homocount)/inputcount
```

## Step 4: Cluster alignments to groups
After filtering alignments, To assemble alignments to DGs and NGs using the crssant.py script, three types of input files are required, `alignfile`, `genesfile` and `bedgraphs`. For more on these parameters, see the explanation below and the bioRxiv preprint referenced at the top of this README. 
```
python crssant.py [-h] [-out OUT] [-cluster CLUSTER] [-n N] [-covlimit COVLIMIT] [-t_o T_O] [-t_eig T_EIG] alignfile genesfile bedgraphs
```

Positional arguments:
* `alignfile`: Path to alignedment file (SAM)
* `genesfile`: Path to gene annotations (BED)
* `bedgraphs`: Path to genome coverage files. Coma separated 2 files for + and - strands

Optional arguments: 
* `-h`: show the help message and exit
* `-out`: path of output folder
* `-cluster`: clustering method, "cliques" or "spectral"(default)
* `-n`: Number of threads. Default is 8
* `-covlimit`: Max coverage to be directly graphed. Default is 1000. 
* `-t_o`: Overlap threshold 0-1 inclusive. Default: 0.5 for "spectral", 0.1 for "cliques"
* `-t_eig`: Eigenratio threshold (positive) for "spectral" only. Default: 5. 


### Required input files
Alignment files `gap1filter.sam` and `trans.sam` are combined to produce `alignfile` sam, keeping one set of header lines at the beginning. The header lines are passed on to output. `genesfile` is a list of all genes in the genome. The start and end for each gene is used to assign the alignments to gene pairs, and determine which alignments correspond to intramolecular structures or intermolecular interactions. The `bedgraphs` files can be produced using the bedtools package, for example: 
```
bedtools genomecov -bg -split -strand + -ibam alignfile_sorted.bam -g chromosome_sizes >alignfile_plus.bedgraph
bedtools genomecov -bg -split -strand - -ibam alignfile_sorted.bam -g chromosome_sizes > alignfile_minus.bedgraph
```

### Clustering methods
The pipeline uses the spectral clustering method to cluster reads into DGs with overlap threshold parameter of `t_o=0.5` and eigenratio threshold of `t_eig=5`. The default spectral clustering method may be operated with different overlap threshold and eigenratio threshold parameters by specifying one or both with the flags `t_o` and `t_eig`, respectively. `t_o` may be any float between 0 and 1, and `t_eig` may be any positive number. Increasing `t_o` tends to result in more DGs that contain fewer reads, and increasing `t_eig` tends to result in fewer DGs containing more reads.

The user may also specify the cliques-finding method for clustering DGs by specifying the clustering flag `cluster` with the `cliques` option. If the cliques-finding method is specified, `t_o` may also be specified, and again may be any float between 0 and 1. The eigenvalue threshold `t_eig` is not used in the cliques-finding method. By default, for the cliques-finding method the overlap threshold is set to 0.1. 

### Output files
After DG clustering, crssant.py verifies that the DGs do not contain any non-overlapping reads, i.e. any reads where the start position of its left arm is greater than or equal to the stop position of the right arm of any other read in the DG. If the DGs do not contain any non-overlapping reads, then the following output files ending in the following are written:

* `_CRSSANT.sam`: SAM file containing alignments that were successfully assigned to DGs, plus DG and NG annotations
* `_CRSSANT_dg.bed`: BED12 file listing all duplex groups. 

The fields of the BED12 file are used according to the [standard definition](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) with the exception of fields 4 and 5. Field 5 `score` is defined as the number of non-continuous alignments in this DG.  Field 4 `name` is defined in the format `gene1,gene2_DGID_covfrac`, where `gene1,gene2` represents the two genes that this DG connects, `DGID` is a numerical ID of the DG, `covfrac` is the confidence of the DG, defined as c / sqrt(a\*b) and
* c = number of reads in a given DG
* a = number of reads overlapping the left arm of the DG
* b = number of reads overlapping the right arm of the DG










## Test

You can test CRSSANT using a collection of Homo sapiens ribosomal RNA (rRNA) test data that we have compiled:

1. Download the compressed folder of [test data](https://github.com/zhipenglu/CRSSANT/tree/master/tests.tar.gz) and decompress using the command `tar -zxvf tests.tar.gz` or by double-clicking on the tar.gz file
2. Specify the path/location where results should be written, e.g. `output`

Run CRSSANT on all rRNA genes in region hs54S:
```
CRSSANT_path/CRSSANT tests/hsrRNA_reads.sam tests/hsrRNA.fa tests/hsrRNA_gene.bed -out output
```
or analyze specific genes, e.g. only reads whose left arms map to gene 5.8S and whose right arms map to gene 28S:
```
CRSSANT_path/CRSSANT tests/hsrRNA_reads.sam tests/hsrRNA.fa tests/hsrRNA_gene.bed -genes 5.8S,28S -out output
```
