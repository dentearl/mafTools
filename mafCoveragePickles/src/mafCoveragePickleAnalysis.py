#!/usr/bin/env python2.7
""" 
mafCoveragePickleGapAnalysis.py
14 November 2011
dent earl

mafCoveragePickleGapAnalysis is a script that operates on a single 
coverage pickle file and extracts information from
the alignments of pairs of sequences. Output is xml which will
contain a 'gaps' tag that contains a comma separated list of all
indels within the pickle for specified pairs. Additionally, 
information on inter-chromosome coverage between pairs is 
provided.

For slightly more detailed examples, check the 
test.mafCoveragePickleGapAnalysis.py unittest.

Comparisons are between a species' chromosome and all other
specified species. So,
for each species S:
    for each chromosome C in S:
        for every species T, T != S:
            paint all positions in C where T aligns (any chrom in T)
"""
##############################
# Copyright (C) 2009-2012 by 
# Dent Earl (dearl@soe.ucsc.edu, dent.earl@gmail.com)
#
#
# ... and other members of the Reconstruction Team of David Haussler's 
# lab (BME Dept. UCSC).
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
##############################
import cPickle
from mafCoveragePickleCreator import readPickle
import numpy
from optparse import OptionParser
import os
from scipy.stats import scoreatpercentile
import re
import sys
import xml.etree.ElementTree as ET

def initOptions(parser):
    parser.add_option('--pickle', dest='pickle',
                      help='input pickle file')
    parser.add_option('--outfile', dest='outfile', default='summary.xml',
                      help='location where outfile will be written. default=%default')
    parser.add_option('--ignore', dest='ignoreList', default='',
                      help='comma separated list of species to ignore. default=%default')
def checkOptions(options, args, parser):
    for k, v in [('pickle', options.pickle), ('outfile', options.outfile),
                 ]:
        if v is None:
            parser.error('specify --%s' % k)
    if not os.path.exists(options.pickle):
        parser.error('--pickle %s does not exist' % options.pickle)
    options.ignoreSet = set(options.ignoreList.split(','))
def maxTheDict(alignmentDict, options):
    maxDict = {}
    for g in alignmentDict:
        maxDict[g] = {}
        for chrom in alignmentDict[g]:
            c = None
            for spp in alignmentDict[g][chrom]:
                ignore = False
                for i in options.ignoreSet:
                    if spp == i:
                        ignore = True
                if ignore:
                    continue
                if c is None:
                    c = alignmentDict[g][chrom][spp]
                else:
                    c = numpy.vstack((c, alignmentDict[g][chrom][spp]))
            maxDict[g][chrom] = {'max': c.max(0)}
    return maxDict
def analyzeAll(alignments, options):
    """ alignments is a multi dict keyed genome1:chrom1:genome2:chrom2: numpy.ndarray()
    where the array is length of genome1:chrom1 and contains 0s in all of the positions
    where genome2:chrom2 does not align.
    """
    cnvsGenomeDict = {}
    for g1 in alignments:
        cnvsGenomeDict[g1] = {}
        for c1 in alignments[g1]:
            for g2 in alignments[g1][c1]:
                analyzeOne(alignments[g1][c1][g2], cnvsGenomeDict[g1], options)
    return cnvsGenomeDict
def analyzeOne(array, cnvs, options):
    """ Take one numpy array and using numpy primitives discover the lengths
    of runs of values and record them in the cnvs dict.
    """
    maxV = numpy.max(array)
    if array.__class__ != numpy.ndarray:
        raise RuntimeError('Type Error, analyzeOne(): array should be numpy.ndarray, not %s\n' 
                           % array.__class__)
    for i in xrange(maxV + 1):
        if i not in cnvs:
            cnvs[i] = numpy.array([], dtype=numpy.uint32)
        cnvs[i] = numpy.concatenate((cnvs[i], findRunsOfValues(array, [i])))
    if 'zerosAndOnes' not in cnvs:
        cnvs['zerosAndOnes'] = numpy.array([], dtype=numpy.uint32)
    cnvs['zerosAndOnes'] = numpy.concatenate((cnvs['zerosAndOnes'], findRunsOfValues(array, [0, 1])))
    if 'greaterThanOne' not in cnvs:
        cnvs['greaterThanOne'] = numpy.array([], dtype=numpy.uint32)
    if maxV > 1:
        if maxV > 2:
            cnvs['greaterThanOne'] = numpy.concatenate((cnvs['greaterThanOne'], 
                                                       findRunsOfValues(array, range(2, maxV + 1))))
        else:
            cnvs['greaterThanOne'] = numpy.concatenate((cnvs['greaterThanOne'], 
                                                       findRunsOfValues(array, [2])))
def findRunsOfValues(array, valuesList):
    """ look through array and find the lengths of runs of values in valuesList.
    Use the values in valuesList to do a logical OR on array and then look for 
    runs of 1.
    """
    if array.__class__ == list:
        array = numpy.array(array, dtype=numpy.uint32)
    if len(valuesList) < 1:
        raise RuntimeError('valuesList must have a length greater than 0, valuesList=%s.' % str(valuesList))
    bits = (array == valuesList[0]) + 0
    for i in xrange(1, len(valuesList)):
        if valuesList[i].__class__ != int:
            raise RuntimeError('Odd, I would have thought valuesList[%d].__class__ == int, not %s' 
                               % (i, valuesList[i].__class__))
        bits = numpy.logical_or(bits, array == valuesList[i]) + 0
    bounded = numpy.hstack(([0], bits, [0])) # puts 0 at start and end of array.
    diffs = numpy.diff(bounded)
    runStarts, = numpy.where(diffs > 0)
    runEnds, = numpy.where(diffs < 0)
    return numpy.array(runEnds - runStarts, dtype=numpy.uint32)
def writeAnalysis(cnvsGenomeDict, combinedDict, alignments, options):
    root = ET.Element('data')
    pc = ET.SubElement(root, 'pairwiseCoverageAnalysis')
    pc.attrib['description'] = ('For every genome g_i, for every chromosome c_j in genome g_i, for every '
                                'other genome g_k (i!=k), report information about the coverage '
                                '(proportion of number of bases with alignments)  of g_k onto c_j.')
    for g1 in alignments:
        g = ET.SubElement(pc, 'targetGenome')
        g.attrib['name'] = g1
        for c1 in alignments[g1]:
            c = ET.SubElement(g, 'targetChromosome')
            c.attrib['name'] = c1
            for g2 in alignments[g1][c1]:
                tlen = len(alignments[g1][c1][g2])
                basesCovered = calcBasesCovered(alignments[g1][c1][g2], options)
                p = ET.SubElement(c, 'partnerGenome')
                p.attrib['name'] = g2
                p.attrib['targetNumberBasesCovered'] = str(basesCovered)
                p.attrib['targetPercentCovered'] = str(basesCovered / float(tlen))
            c.attrib['targetLength'] = str(tlen)

    e = ET.SubElement(root, 'copyNumberAnalysis')
    f = ET.SubElement(e, 'eachArrayWalked')
    f.attrib['description'] = ('For each genome, for each chromosome, each alignment array to '
                               'another genome is individually walked and the lengths of runs '
                               'of values are recorded in a dictionary keyed on length of run '
                               'and valued on number of times observed. ')
    for g1 in cnvsGenomeDict:
        g = ET.SubElement(f, 'genome')
        g.attrib['name'] = g1
        sortedCnvs = sorted(cnvsGenomeDict[g1], key=lambda x: x)
        for s in sortedCnvs:
            c = ET.SubElement(g, 'copy_%s' % str(s))
            if len(cnvsGenomeDict[g1][s]) == 0:
                c.attrib['max'] = str(float('nan'))
                c.attrib['lastQuart'] = str(float('nan'))
                c.attrib['median'] = str(float('nan'))
                c.attrib['mean'] = str(float('nan'))
                c.attrib['firstQuart'] = str(float('nan'))
                c.attrib['min'] = str(float('nan'))
                c.attrib['n'] = str(0)
                c.attrib['sum'] = str(0)
            else:
                c.attrib['max'] = str(max(cnvsGenomeDict[g1][s]))
                c.attrib['lastQuart'] = str(scoreatpercentile(cnvsGenomeDict[g1][s], 75))
                c.attrib['median'] = str(numpy.median(cnvsGenomeDict[g1][s]))
                c.attrib['mean'] = str(numpy.mean(cnvsGenomeDict[g1][s]))
                c.attrib['firstQuart'] = str(scoreatpercentile(cnvsGenomeDict[g1][s], 25))
                c.attrib['min'] = str(min(cnvsGenomeDict[g1][s]))
                c.attrib['n'] = str(len(cnvsGenomeDict[g1][s]))
                c.attrib['sum'] = str(sum(cnvsGenomeDict[g1][s]))
            cnvsGenomeDict[g1][s].sort()
            c.text = ','.join(map(str, cnvsGenomeDict[g1][s]))
    f = ET.SubElement(e, 'arraysCombinedIntoOneThenWalked')
    f.attrib['description'] = ('For each genome, for each chromosome, all alignment arrays are combined '
                               'into a single array by way of a max performed on each '
                               'column. This array is then walked looking for runs of values exactly as '
                               'in the "eachArrayWalked" analysis.')
    for g1 in cnvsGenomeDict:
        g = ET.SubElement(f, 'genome')
        g.attrib['name'] = g1
        sortedCnvs = sorted(combinedDict[g1], key=lambda x: x)
        for s in sortedCnvs:
            c = ET.SubElement(g, 'copy_%s' % str(s))
            if len(combinedDict[g1][s]) == 0:
                c.attrib['max'] = str(float('nan'))
                c.attrib['lastQuart'] = str(float('nan'))
                c.attrib['median'] = str(float('nan'))
                c.attrib['mean'] = str(float('nan'))
                c.attrib['firstQuart'] = str(float('nan'))
                c.attrib['min'] = str(float('nan'))
                c.attrib['n'] = str(0)
                c.attrib['sum'] = str(0)
            else:
                c.attrib['max'] = str(max(combinedDict[g1][s]))
                c.attrib['lastQuart'] = str(scoreatpercentile(combinedDict[g1][s], 75))
                c.attrib['median'] = str(numpy.median(combinedDict[g1][s]))
                c.attrib['mean'] = str(numpy.mean(combinedDict[g1][s]))
                c.attrib['firstQuart'] = str(scoreatpercentile(combinedDict[g1][s], 25))
                c.attrib['min'] = str(min(combinedDict[g1][s]))
                c.attrib['n'] = str(len(combinedDict[g1][s]))
                c.attrib['sum'] = str(sum(combinedDict[g1][s]))
            combinedDict[g1][s].sort()
            c.text = ','.join(map(str, combinedDict[g1][s]))
    info = ET.ElementTree(root)
    info.write(options.outfile)
def calcBasesCovered(array, options):
    """ Function-ized for unittesting.
    """
    return len(numpy.nonzero(array)[0])
def main():
    usage = ('usage: %prog --pickle path/to/file.pickle\n\n'
             '%prog is a script that operates on a single maf\n'
             'mafCoveragePickleGapAnalysis is a script that operates on a single \n'
             'coverage pickle file and extracts information from\n'
             'the alignments of pairs of sequences. Output is xml which will\n'
             'contain a \'gaps\' tag that contains a comma separated list of all\n'
             'indels within the pickle for specified pairs. Additionally, \n'
             'information on inter-chromosome coverage between pairs is \n'
             'provided.\n\n'
             'For slightly more detailed examples, check the \n'
             'test.mafCoveragePickleGapAnalysis.py unittest.\n\n'
             'Comparisons are between a species\' chromosome and all other\n'
             'specified species. So,\n'
             'for each species S:\n'
             '    for each chromosome C in S:\n'
             '        for every species T, T != S:\n'
             '            paint all positions in C where T aligns (any chrom in T)'
             )
    parser = OptionParser(usage = usage)
    initOptions(parser)
    options, args = parser.parse_args()
    checkOptions(options, args, parser)
    
    alignmentDict = readPickle(options.pickle, options)
    cnvsGenomeDict = analyzeAll(alignmentDict, options)
    maxedDict = maxTheDict(alignmentDict, options)
    combinedDict = analyzeAll(maxedDict, options)
    
    writeAnalysis(cnvsGenomeDict, combinedDict, alignmentDict, options)
    
if __name__ == '__main__':
    main()
