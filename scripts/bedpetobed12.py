"""
bedpetobed12.py. Zhipeng Lu. 2020, zhipengluchina@gmail.com

This script converts "convertable" bedpe records to bed12 format, which means
in the record in bedpe satisfies (chrom1 == chrom2) and (strand1 == strand2). 

bed12 (12 fields): 
chrom chromStart chromEnd name score strand thickStart thickEnd
itemRgb blockCount blockSizes blockStarts

bedpe (10 fields): 
chrom1 chromStart1 chromEnd1 chrom2 chromStart2 chromEnd2
name score strand1 strand2
"""

import sys

if len(sys.argv) < 3:
    print("Usage: python bedpetobed12.py bedpefile bed12file")
    sys.exit()

bedpefile = open(sys.argv[1], 'r')
bed12file = open(sys.argv[2], 'w')
outlist = []

for line in bedpefile:
    record = line.strip('\n').split()
    chrom = record[0]
    chromStart1, chromEnd1 = int(record[1]), int(record[2])
    chromStart2, chromEnd2 = int(record[4]), int(record[5])
    name, score, strand = record[6], record[7], record[8]
    itemRgb = "0,0,0"
    blockSizes = str(chromEnd1-chromStart1)+","+str(chromEnd2-chromStart2)
    blockStarts = "0," + str(chromStart2-chromStart1)
    outrecord = [chrom, str(chromStart1), str(chromEnd2), name, score, \
                 strand, str(chromStart1), str(chromStart1), itemRgb, "2", \
                 blockSizes, blockStarts]
    outlist.append("\t".join(outrecord) + "\n")

bed12file.write(''.join(outlist))
bedpefile.close()
bed12file.close()
