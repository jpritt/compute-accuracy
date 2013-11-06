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
def read_data(filename, form, junctions, wgt):
    print 'Reading ' + filename

    with open(filename, 'r') as tsv:
        rows = []
        for line in tsv:
            row = line.strip().split('\t')
            if (form == 'sam'):
                chr_name = row[2]
                if chr_name in chromosomes:
                    start = int(row[3])

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

                    #print origPattern + ' -> ' + str(junctionOffsets)
                    
                    for off in junctionOffsets:
                        i = off + start - 1
                        junctions[i] = 1
                #else:
                #    print 'Unknown chromosome ' + row[2]
            elif (form == 'bed'):
                chr_name = row[0]
                if chr_name == 'dmel_mitochondrion_genome':
                    chr_name = 'M'
                if chr_name in chromosomes:
                    start = int(row[1])
                    junctionLens = [int(s) for s in string.split(row[10], ',')]
                    junctionOffsets = [int(s) for s in string.split(row[11], ',')]

                    for i in xrange(len(junctionOffsets)):
                        off = junctionOffsets[i]
                        length = junctionLens[i]
                        if i > 0:
                            junctions[off+start] = 1
                        if i < (len(junctionOffsets)-1):
                            junctions[off+start+length] = 1
                #else:
                #    print 'Unknown chromosome ' + chr_name

            else:
                chr_name = row[0]
                if chr_name == 'dmel_mitochondrion_genome':
                    chr_name = 'M'
 
                if row[2] == 'exon':
                    junctions[int(row[3])-1] = 1
                    junctions[int(row[4])] = 1


    return junctions

def findJunctions(filename, radius, threshold):
    total_length = 0
    for [a,b] in chr_lengths:
        total_length += b
    cov = np.zeros(total_length)

    with open(filename, 'r') as tsv:
        for line in tsv:
            row = line.strip().split('\t')

            chr_name = row[2]
            if chr_name in chromosomes:
                start = int(row[3])
                
                pattern = row[5]
                sectionLens = []
                sectionOffsets = []
                currOffset = 0

                match=  re.search("\D", pattern)
                while match:
                    index = match.start()
                    if pattern[index] == 'M':
                        sectionOffsets.append(currOffset)
                        sectionLens.append(int(pattern[:index]))
                    currOffset += int(''.join(pattern[:index]))
                    pattern = pattern[index+1:]
                    match = re.search("\D", pattern)
                for i in xrange(len(sectionLens)):
                    for j in xrange(sectionLens[i]):
                        index = getIndex(start+sectionOffsets[i]+j, chr_name)
                        cov[index] += 1

    # take derivative at each base
    #wgt = float(1) / (2*radius+1)
    #for i in xrange(total_length-1):
    #    for j in xrange(-radius,radius+1):
    #        if i+j >= 0 and i+j < total_length:
    #            cov[i+j] += wgt*(cov[i+1] - cov[i])

    for i in xrange(total_length-1):
        cov[i] = cov[i+1] - cov[i]

    junctions = dict()

    for i in xrange(1, total_length-1):
        v = cov[i]
        if abs(v) > threshold:
            if cov[i] > cov[i-1] and cov[i] > cov[i+1] and cov[i] > 0:
                junctions[i] = 1
            elif cov[i] < cov[i-1] and cov[i] < cov[i+1] and cov[i] < 0:
                junctions[i] = 1
        
    return junctions


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


if len(sys.argv) != 3:
    print "Usage: ./comp_coverage.py actual predicted"
    print "  .sam .bed and .gtf files are supported"

actual = sys.argv[1]
predicted = sys.argv[2]

format1 = actual[-3:len(actual)]
format2 = predicted[-3:len(predicted)]

if (format1 != 'sam' and format1 != 'bed' and format1 != 'gtf') or (format2 != 'sam' and format2 != 'bed' and format2 != 'gtf'):
    print 'Only .sam, .bed, and .gtf files are supported'


#print 'Reading actual alignments'
j_actual = dict()
j_actual = read_data(actual, format1, j_actual, 1)

print len(j_actual)


for thresh in xrange(11):
    print 'Calculating predicted alignments for threshold ' + str(thresh)
    j_predicted = findJunctions(predicted, 2, thresh)

    print '  ' + len(j_predicted)
    tp = 0
    fp = 0
    fn = 0
    tn = 0

    jarray_actual = np.array([])
    jarray_predicted = np.array([])

    for i in j_predicted:
        jarray_predicted = np.append(jarray_predicted, j_predicted[i])
        if i in j_actual:
            tp += 1
            jarray_actual = np.append(jarray_actual, j_actual[i])
        else:
            fp += 1
            jarray_actual = np.append(jarray_actual, 0)
    for i in j_actual:
        if not i in j_predicted:
            fn += 1
            jarray_actual = np.append(jarray_actual, j_actual[i])
            jarray_predicted = np.append(jarray_predicted, 0)

    for (name, length) in chr_lengths:
        tn += length
    tn = tn - tp - fp - fn

    total = len(jarray_actual)

    print '  True positives:\t' + str(tp)
    print '  False positives:\t' + str(fp) 
    print '  False negatives:\t' + str(fn)
    #print 'True negatives:\t\t' + str(tn)
    #print '\n'

    #print 'Sensitivity:\t\t' + str(float(tp) / (tp+fn))
    #print 'Specificity:\t\t' + str(float(tn) / (tn+fp))
    #print 'Accuracy:\t\t' + str(float(tp+tn) / (tp+fp+fn+tn))
    print '  Precision:\t\t' + str(float(tp) / (tp+fp)) 
    print '  Recall:\t\t' + str(float(tp) / (tp+fn))
    #print '\n'

    dist = np.linalg.norm(np.subtract(jarray_actual, jarray_predicted))

    len_actual = np.linalg.norm(jarray_actual)
    len_predicted = np.linalg.norm(jarray_predicted)
    jarray_actual = jarray_actual / len_actual
    jarray_predicted = jarray_predicted / len_predicted

    correlation = np.correlate(jarray_actual, jarray_predicted)

    #print '  Distance:\t' + str(dist)
    print '  Correlation:\t' + str(correlation)


# plot all points
#print 'Plotting'
#matplotlib.pyplot.scatter(jarray_actual, jarray_predicted)
#matplotlib.pyplot.show()
print 'Done!'
