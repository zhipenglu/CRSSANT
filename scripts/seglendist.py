"""
seglendist.py. Zhipeng Lu, 2021, zhipengluchina@gmail.com
This script calculates the size distribution of each segment in STAR alignments. 

Requires python3.

Example command for creating the list:
python seglendist.py align.sam sam align_seglen.list
Example command for making the distribution figure: 
python seglendist.py align_seglen.list list align_seglen.pdf
"""

import matplotlib.pyplot as plt
import sys, re, matplotlib, matplotlib.pyplot as plt

if len(sys.argv)< 4:
    print("Usage: python armsizedist.py inputfile filetype outputfile")
    print("filetype: sam or list")
    print("when input is sam, output is list")
    print("when input is list, output is the pdf figure")
    sys.exit()

inputfile = open(sys.argv[1], 'r')
filetype = sys.argv[2]


##part1: save the size list as space delimited numbers in one line in a file. 
if filetype == "sam":
    outputfile = open(sys.argv[3], 'w')
    sizelist = []
    for line in inputfile:
        if line[0] == "@": continue
        CIGAR = line.split()[5]
        segs=[i.rstrip('0123456789') for i in CIGAR.split('N')]
        Mlens=[sum([int(i[:-1])for i in re.findall('\d+[M=X]',s)])for s in segs]
        sizelist += Mlens
    #save the size list as space delimited numbers in one line in a file. 
    sizestring = ' '.join([str(size) for size in sizelist])
    outputfile.write(sizestring)
    outputfile.close()
    inputfile.close()

##part2: plot the arm size distribution
elif filetype == "list":
    sizelist = [int(num) for num in inputfile.readline().split()]
    fig, ax = plt.subplots(figsize=(2.5,1.5))
    plt.subplots_adjust(bottom=0.32)
    plt.subplots_adjust(top=0.95)
    plt.subplots_adjust(left=0.25)
    plt.subplots_adjust(right=0.95)

    n, bins, patches = plt.hist(sizelist, [x+0.5 for x in range(0,100)], \
                                histtype='step', cumulative=True, normed=1) 
    plt.setp(patches, 'facecolor', 'b', 'alpha', 0.75)
    plt.xlim(-10, 50)
    plt.ylim(0, 1.1)
    max_yticks = 5
    yloc = plt.MaxNLocator(max_yticks)
    ax.yaxis.set_major_locator(yloc)
    ax.set_xlabel("segment length (nt)") #, fontsize=15
    ax.set_ylabel("cumulative frequency") #, fontsize=15
    plt.savefig(sys.argv[3])
    plt.show()
           

else:
    print("Input file type is wrong. Please use sam or list")
    sys.exit()
