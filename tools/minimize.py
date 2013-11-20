#! /usr/bin/env python
import sys
import csv
import math
import re
import string
import numpy as np
import pickle

# Convert tophat output to a smaller data structure
# The minimized structure consists of the coverage vector for all reads and a shorter representation of spanning reads

# Read the alignment data from a sam file
# Update the coverage arrays
def minimize(filename, form):
    print 'Reading ' + filename

    spanning = []
    bases = sum([y for [x,y] in chr_lengths])
    print bases
    coverage = np.zeros(bases)

    with open(filename, 'r') as tsv:
        rows = []
        for line in tsv:
            row = line.strip().split('\t')
            if (form == 'sam'):
                chr_name = row[2]
                if chr_name in chromosomes:
                    start = getIndex(int(row[3]), chr_name)

                    pattern = row[5]
                    junctionOffsets = []
                    currOffset = 0

                    origPattern = pattern

                    match = re.search("\D", pattern)
                    while match and match.start() < len(pattern)-1:
                        index = match.start()
                        currOffset += int(''.join(pattern[:index]))
                        pattern = pattern[index+1:]
                        junctionOffsets.append(currOffset)
                        match = re.search("\D", pattern)

                    end = currOffset + int(''.join(pattern[:match.start()]))

                    introns = []
                    if len(junctionOffsets) > 0:
                        for i in xrange(0, len(junctionOffsets), 2):
                            introns.append((junctionOffsets[i] + start, junctionOffsets[i+1] + start))
                        d0 = junctionOffsets[0]
                        df = end - junctionOffsets[-1]
                        spanning.append([d0, df, introns])
                    else:
                        for i in xrange(start,start+end):
                            coverage[i] += 1
    return coverage, spanning


chromosomes = ['2L', '2R', '3L', '3R', '4', 'M', 'X', '2LHet', '2RHet', '3LHet', '3RHet', 'XHet', 'YHet', 'U', 'Uextra']
chr_lengths = [['2L', 23011544],
               ['2R', 21146708],
               ['3L', 24543557],
               ['3R', 27905053],
               ['4', 1351857],
               ['M', 19517],
               ['X', 22422827],
               ['2LHet', 368872],
               ['2RHet', 3288761],
               ['3LHet', 2555491],
               ['3RHet', 2517507],
               ['XHet', 204112],
               ['YHet', 347038],
               ['U', 10049037],
               ['Uextra', 29004656]]

#chromosomes = ['2L']
#chr_lengths = [['2L', 23011544]]

# Return the position in the 1d arry of the nt at the given position in a chromosome
def getIndex(position, chromosome):
    index = 0
    for [a,b] in chr_lengths:
        if a == chromosome:
            return index + position
        else:
            index += b


# Initialize empty coverage arrays
j_actual = dict()
j_predicted = dict()

if len(sys.argv) != 2:
    print "Usage: ./minimize.py alignments.sam"

alignments = sys.argv[1]

form = alignments[-3:len(alignments)]

if (format != 'sam'):
    print 'Only .sam files are supported'

coverage, spanning = minimize(alignments, form)

print 'Writing coverage file'
pickle.dump(' '.join(coverage), 'coverage.txt')

print 'Writing spanning file'
f = open('spanning.txt', 'w')
for row in spanning:
    f.write(str(spanning[0]) + ' ' + str(spanning[1]))
    for (x,y) in spanning[2]:
        f.write(' ' + str(x) + ' ' + str(y))
    f.write('\n')
f.close()
