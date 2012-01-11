#!/usr/bin/env python2.6
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
import mafCoveragePicklePlotter as mcpp
import numpy
from optparse import OptionParser
import os
import re
import xml.etree.ElementTree as ET

def initOptions(parser):
    parser.add_option('--pickle', dest = 'pickle',
                      help = 'input pickle file')
    parser.add_option('--outfile', dest = 'outfile', default = 'summary.xml',
                      help = 'location where outfile will be written default = %default')
    parser.add_option('--noEdges', dest = 'noEdges', action = 'store_true', default = False,
                      help = ('Gaps redefined to occur only between areas of alignments. Has '
                              'the effect of turning off gaps that occur at the edges of chromosomes. '
                              'default = %default.'))

def checkOptions(options, args, parser):
    for k, v in [('pickle', options.pickle), ('outfile', options.outfile),
                 ]:
        if v is None:
            parser.error('specify --%s' % k)
    if not os.path.exists(options.pickle):
        parser.error('--pickle %s does not exist' % options.pickle)

def analyzeAll(alignments, options):
    """ alignments is a multi dict keyed genome1:chrom1:genome2:chrom2: numpy.array()
    where the array is length of genome1:chrom1 and contains 0s in all of the positions
    where genome2:chrom2 does not align.
    """
    gaps = []
    for g1 in alignments:
        for c1 in alignments[g1]:
            for g2 in alignments[g1][c1]:
                gaps += analyzeOne(alignments[g1][c1][g2], options)
    return gaps

def analyzeOne(array, options):
    """ array is a numpy.array(), 1 by n. 0s indicate no alignment,
    x > 0 indicate that x bases aligned to that position.
    """
    result = []
    count = 0
    startEdge = True
    for a in array:
        if a == 0:
            count += 1
        else:
            if options.noEdges and startEdge:
                startEdge = False
                count = 0
                continue
            if count != 0:
                result.append(count)
                count = 0
    if not options.noEdges:
        if count != 0:
            result.append(count)
    return result

def writeAnalysis(gapsList, alignments, options):
    root = ET.Element('data')
    pc = ET.SubElement(root, 'pairwiseCoverage')
    for g1 in alignments:
        for c1 in alignments[g1]:
            for g2 in alignments[g1][c1]:
                tlen = len(alignments[g1][c1][g2])
                basesCovered = calcBasesCovered(alignments[g1][c1][g2], options)
                e = ET.SubElement(pc, 'coverage')
                e.attrib['targetGenome'] = g1
                e.attrib['targetChromosome'] = c1
                e.attrib['queryGenome'] = g2
                e.attrib['targetLength'] = str(tlen)
                e.attrib['targetNumberBasesCovered'] = str(basesCovered)
                e.attrib['targetPercentCovered'] = str(basesCovered / float(tlen))
    e = ET.SubElement(root, 'gaps')
    e.text = ','.join(map(str, gapsList))
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
    
    alignmentDict = mcpp.readPickle(options.pickle, options)
    gapsList = analyzeAll(alignmentDict, options)
    writeAnalysis(gapsList, alignmentDict, options)
    
if __name__ == '__main__':
    main()
