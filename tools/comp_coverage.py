#! /usr/bin/env python
import sys
import csv
import math
import re
import string
import numpy as np
import matplotlib
import pylab

# Description goes here

# Read the alignment data from a sam file
# Update the coverage arrays
def read_sam(filename, form, cov):
    print 'Reading ' + filename

    hits = 0
    with open(filename, 'r') as tsv:
        rows = []
        for line in tsv:
            row = line.strip().split('\t')
            wgt = 1
            if (form == 'sam'):
                chr_name = row[2]
                if chr_name in chromosomes:
                    start = int(row[3])

                    pattern = row[5]
                    sectionLens = []
                    sectionOffsets = []
                    currOffset = 0

                    match = re.search("\D", pattern)
                    while match:
                        index = match.start()
                        if pattern[index] == 'M':
                            sectionOffsets.append(currOffset)
                            sectionLens.append(int(pattern[:index]))
                        currOffset += int(''.join(pattern[:index]))
                        pattern = pattern[index+1:]
                        match = re.search("\D", pattern)

                    for i in row[11:len(row)]:
                        if i[0:5] == 'NH:i:':
                            wgt = 1 / int(i[5:len(i)])
            else:
                chr_name = row[0]
                if chr_name == 'dmel_mitochondrion_genome':
                    chr_name = 'M'
                if chr_name in chromosomes:
                    start = int(row[1])
                    
                    sectionLens = [int(s) for s in string.split(row[10], ',')]
                    sectionOffsets = [int(s) for s in string.split(row[11], ',')]
                #else:
                #    print 'Unknown chromosome ' + chr_name

            if chr_name in chromosomes:
                for i in xrange(len(sectionLens)):
                    for j in xrange(sectionLens[i]):
                        index = getIndex(start+sectionOffsets[i]+j, chr_name)
                        cov[index] += wgt
                    hits += sectionLens[i];

    return cov,hits;

    

#chromosomes = ['2L', '2R', '3L', '3R', '4', 'M', 'X', '2LHet', '2RHet', '3LHet', '3RHet', 'XHet', 'YHet', 'U', 'Uextra']
#chr_lengths = [['2L', 23011544],
#               ['2R', 21146708],
#               ['3L', 24543557],
#               ['3R', 27905053],
#               ['4', 1351857],
#               ['M', 19517],
#               ['X', 22422827],
#               ['2LHet', 368872],
#               ['2RHet', 3288761],
#               ['3LHet', 2555491],
#               ['3RHet', 2517507],
#               ['XHet', 204112],
#               ['YHet', 347038],
#               ['U', 10049037],
#               ['Uextra', 29004656]]

chromosomes = ['2L']
chr_lengths = [['2L', 23011544]]

# Return the position in the 1d arry of the nt at the given position in a chromosome
def getIndex(position, chromosome):
    index = 0
    for [a,b] in chr_lengths:
        if a == chromosome:
            return index + position
        else:
            index += b


# Construct coverage arrays
print 'Initializing coverage arrays'
total_length = 0
for [a, b] in chr_lengths:
    total_length += b

coverage_actual = np.zeros(total_length)
coverage_predicted = np.zeros(total_length)

if len(sys.argv) != 3:
    print "Usage: ./comp_coverage.py actual predicted"
    print "  .sam and .bed files are supported"

actual = sys.argv[1]
predicted = sys.argv[2]

format1 = actual[-3:len(actual)]
format2 = predicted[-3:len(predicted)]

if (format1 != 'sam' and format1 != 'bed') or (format2 != 'sam' and format2 != 'bed'):
    print 'Only .sam and .bed files are supported'


print 'Reading actual alignments'
coverage_actual, hits_actual = read_sam(actual, format1, coverage_actual)
print str(hits_actual) + ' actual hits'

print 'Calculating length'
len_actual = np.linalg.norm(coverage_actual)
print 'Actual Length: ' + str(len_actual)

print 'Reading predicted alignments'
coverage_predicted, hits_predicted = read_sam(predicted, format2, coverage_predicted)
print str(hits_predicted) + ' predicted hits'

print 'Calculating error'
len_predicted = np.linalg.norm(coverage_predicted)
print 'Predicted Length: ' + str(len_predicted)

# normalize
coverage_actual = coverage_actual / len_actual
coverage_predicted = coverage_predicted / len_predicted

dist = np.linalg.norm(np.subtract(coverage_actual, coverage_predicted))
print 'Distance: ' + str(dist)

correlation = np.correlate(coverage_actual, coverage_predicted)

print 'Distance:\t' + str(dist)
print 'Correlation:\t' + str(correlation)


# plot every 10th point, ignoring all (0,0) points
print 'Plotting'
xs = coverage_actual[0:len(coverage_actual):10]
ys = coverage_predicted[0:len(coverage_predicted):10]
matplotlib.pyplot.scatter(xs+ys, xs-ys, alpha=0.5)
matplotlib.pyplot.show()
print 'Done!'
