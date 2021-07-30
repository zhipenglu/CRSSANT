"""
Zhipeng Lu, 2020, zhipengluchina@gmail.com
Requires python 3. 

gaplendist.py, a simple script to plot the length of gaps (N). The input file
can be any type of SAM. For example, from the direct output of STAR
we only have the distribution from normally gapped reads. 

Command for creating the list:
python gaplendist.py input.sam sam input_gaplen.list
Command for creating the distribution pdf file:
python gaplendist.py input_gaplen.list list input_gaplen.pdf
"""

import sys, re, matplotlib
import matplotlib.pyplot as plt

if len(sys.argv)< 5:
    print("Usage: python gaplendist.py inputfile filetype outfile gap")
    print("type: 'sam' file or 'list' of lengths")
    print("when input is sam, output is list")
    print("when input if list, output is pdf figure")
    print("gap: all or min (smallest one only for each alignment)")
    sys.exit()
inputfile = open(sys.argv[1], 'r')
filetype = sys.argv[2]
gap = sys.argv[4]


##part1: save the size list as space delimited numbers in one line in a file. 
if filetype == "sam":
    outputfile = open(sys.argv[3], 'w')
    sizelist = [] #a file of integer numbers
    for line in inputfile:
        if line[0] == "@": continue
        align = line.split()
        cigar = align[5]
        if gap=="min":sizelist+=[min([i[:-1]for i in re.findall('\d+N',cigar)])]
        if gap=="all":sizelist+=[i[:-1]for i in re.findall('\d+N',cigar)]
    print "Total number of gaps", len(sizelist)
    sizestring = ' '.join([str(size) for size in sizelist])
    outputfile.write(sizestring)
    outputfile.close()
    inputfile.close()


##part2: plot the gap size distribution
elif filetype == "list":
    sizelist = [int(i) for i in inputfile.readline().split()]
    fig, ax = plt.subplots(figsize=(2.5,1.5))
    plt.subplots_adjust(bottom=0.32)
    plt.subplots_adjust(top=0.95)
    plt.subplots_adjust(left=0.25)
    plt.subplots_adjust(right=0.95)
    n, bins, patches = plt.hist(sizelist, [x+0.5 for x in range(0,100)], \
                                histtype='step', cumulative=True, normed=1)
    plt.setp(patches, 'facecolor', 'b', 'alpha', 0.75)
    plt.xlim(-10, 100)
    plt.ylim(0, 1.1)

    max_yticks = 5
    yloc = plt.MaxNLocator(max_yticks)
    ax.yaxis.set_major_locator(yloc)
    ax.set_xlabel("gap length (nt)") #, fontsize=15
    ax.set_ylabel("frequency") #, fontsize=15
    plt.savefig(sys.argv[3])
    plt.show()

else:
    print("Input file type is wrong. Please use sam or list")
    sys.exit()
