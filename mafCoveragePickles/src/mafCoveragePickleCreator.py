#!/usr/bin/env python2.6
""" 
mafCoveragePickleCreator.py
10 January 2011
dent earl

mafCoveragePickleCreator is a script that operates on a single maf
(multiple alignment format) file and extracts information from
the alignments of pairs of sequences. Output is a pickle.
The maf sequence name field must have the format
species.chromosome
e.g.
hg19.chr22

For slightly more detailed examples, check the 
test.mafCoveragePickleCreator.py unittest.

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
import numpy
from optparse import OptionParser
import os
import re

seqRegex = r'[ACGTUMRWSYKVHDBXN]+'
seqPat = re.compile(seqRegex)

class MafLine:
    def __init__(self, genome = '', chrom = '', start = -1, 
                 seqLength = -1, totalLength = -1, strand = 0,
                 sequence = -1, lineno = -1):
        self.genome = genome
        self.chrom = chrom
        self.start = start
        self.seqLength = seqLength
        self.strand = strand
        self.totalLength = totalLength
        self.sequence = sequence
        self.lineno = lineno

def initOptions(parser):
    parser.add_option('--maf', dest = 'maf', 
                      help = 'input maf file')
    parser.add_option('--species', dest = 'species', 
                      help = 'comma separated list of species names to include in output')
    parser.add_option('--pickle', dest = 'pickle',
                      help = ('location where python pickle will be written, for use '
                      'with downstream analyses. By default this file is not created.'))

def checkOptions(options, args, parser):
    for k, v in [('maf', options.maf), ('pickle', options.pickle),
                 ('species', options.species)]:
        if v is None:
            parser.error('specify --%s' % k)
    if not os.path.exists(options.maf):
        parser.error('--maf %s does not exist' % options.maf)
    if options.species is None:
        options.speciesList = set()
    else:
        options.speciesList = set(options.species.split(','))

def readMaf(filename, options):
    """ read a given maf file and create an array 'alignments'
    which is a multilayered dict:
    alignments[a.genome][a.chrom][b.genome][b.chrom] = numpy.zeros(a.totalLength)
    where the numpy array represents the locations where genome b maps onto genome a
    """
    alignments = {}
    f = open(filename, 'r')
    
    blocks = []
    mafLineList = []
    
    namepat = re.compile(r'(.+?)\.(.*)')
    
    for lineno, line in enumerate(f, 1):
        line = line.strip()
        if line.startswith('s'):
            ml = extractMafLine(namepat, line, lineno, options)
            if ml is not None:
                mafLineList.append(ml)
        else:
            if len(mafLineList) > 0:
                addBlockPairs(alignments, mafLineList, options)
            mafLineList = []
    
    if len(mafLineList) > 0:
        addBlockPairs(alignments, mafLineList, options)
        mafLineList = []
    f.close()
    return alignments

def extractMafLine(namepat, line, lineno, options):
    """ parse a given line from a maf file into a 
    MafLine object.
    """
    data = line.split()
    if len(data) != 7:
        raise RuntimeError('maf line with incorrect number of fields on lineno %d: %s' % (lineno, line))
    ml = MafLine()
    m = re.match(namepat, data[1])
    if m is None:
        raise RuntimeError('maf sequence line where name field has no chr on lineno %d: %s' % (lineno, line))
    ml.genome = m.group(1)
    if ml.genome not in options.speciesList:
        return None
    ml.chrom = m.group(2)
    ml.start  = int(data[2])
    ml.seqLength = int(data[3])
    ml.strand = int(data[4] + '1')
    ml.totalLength = int(data[5])
    ml.sequence = data[6]
    ml.lineno = lineno
    return ml

def addBlockPairs(alignments, mafLineList, options):
    """ loop through all pairs of mafLines in the block list, store
    the discovered alignments in the alignments dict.
    """
    for a in mafLineList:
        for b in mafLineList:
            if a.genome == b.genome:
                # we skip self-alignments
                continue
            if a.genome not in alignments:
                alignments[a.genome] = {}
            if a.chrom not in alignments[a.genome]:
                alignments[a.genome][a.chrom] = {}
            if b.genome not in alignments[a.genome][a.chrom]:
                # each 'a' chrom will have an array for each 'b' genome, so long
                # as they appear together in an alignment block.
                alignments[a.genome][a.chrom][b.genome] = numpy.zeros(a.totalLength, 
                                                                      dtype = numpy.uint16)
            addBlocksToArray(alignments[a.genome][a.chrom][b.genome], a, b)
    
    # explicitly throw this away to help with memory
    mafLineList = []

def addBlocksToArray(array, a, b):
    """ array is a numpy array of columns in genome/chrom a that 
    are aligned to any position in genome/chrom b.
    i.e. the length of 'genome/chrom a' = length of 'array'
    a and b are both mafline objects.
    """
    global seqPat
    
    aList = [] # ranges as tuples
    bList = [] 
    # offsets allow us to keep track of the 
    # within genome coordinates of the sequence field.
    # they are offsets from the sequence start field integer
    # and are in the sequence strand direction
    aOffsets = [] 
    bOffsets = [] 
    
    ##########
    # find all of the contiguous intervals of gapless sequence
    # in the sequence fields of a and b.
    # ranges are mathematically [start, end), i.e. half-open.
    for m in re.finditer(seqPat, a.sequence):
        aList.append((m.start(), m.end()))
        if len(aOffsets) == 0:
            aOffsets.append(0)
        else:
            aOffsets.append(aOffsets[-1] + m.end() - m.start())
    for m in re.finditer(seqPat, b.sequence):
        bList.append((m.start(), m.end()))
        if len(bOffsets) == 0:
            bOffsets.append(0)
        else:
            bOffsets.append(bOffsets[-1] + m.end() - m.start())
    
    i, j = 0, 0 # index
    prevI, prevJ = -1, -1
    while i < len(aList) and j < len(bList):
        if (i, j) == (prevI, prevJ):
            print locals()
            raise RuntimeError('Infinite loop detected')
        prevI, prevJ = i, j
        if aList[i][1] <= bList[j][0]:
            # i ---
            # j     ---
            # no overlap
            i += 1
            continue
        if bList[j][1] <= aList[i][0]:
            # i     ---
            # j ---
            # no overlap
            j += 1
            continue
        if aList[i][1] < bList[j][1] and aList[i][1] > bList[j][0] and aList[i][0] < bList[j][0]:
            # i ----
            # j   ----
            # overlap is bList[j][0]..aList[i][1]
            # forgive the overly verbose indexing, it makes it easier to understand.
            forStart = a.start + aOffsets[i] - aList[i][0] + bList[j][0]
            forEnd   = a.start + aOffsets[i] - aList[i][0] + aList[i][1]
            incrementArrayIndices(array, a.totalLength, forStart, forEnd, a.strand)
            i += 1
            continue
        if bList[j][1] < aList[i][1] and bList[j][1] > aList[i][0] and bList[j][0] < aList[i][0]:
            # i   ----
            # j ----
            # overlap is aList[i][0]..bList[j][1]
            forStart = a.start + aOffsets[i] - aList[i][0] + aList[i][0]
            forEnd   = a.start + aOffsets[i] - aList[i][0] + bList[j][1]
            incrementArrayIndices(array, a.totalLength, forStart, forEnd, a.strand)
            j += 1
            continue
        if aList[i][0] >= bList[j][0] and aList[i][1] <= bList[j][1]:
            # i   ----
            # j --------
            # overlap is aList[i][0]..aList[i][1]
            forStart = a.start + aOffsets[i] - aList[i][0] + aList[i][0]
            forEnd   = a.start + aOffsets[i] - aList[i][0] + aList[i][1]
            incrementArrayIndices(array, a.totalLength, forStart, forEnd, a.strand)
            i += 1
            continue
        if bList[j][0] >= aList[i][0] and bList[j][1] <= aList[i][1]:
            # i --------
            # j   ----
            forStart = a.start + aOffsets[i] - aList[i][0] + bList[j][0]
            forEnd   = a.start + aOffsets[i] - aList[i][0] + bList[j][1]
            incrementArrayIndices(array, a.totalLength, forStart, forEnd, a.strand)
            j += 1
            continue
        print locals()
        raise RuntimeError('Unanticipated condition.')

def incrementArrayIndices(array, totalLength, forStart, forEnd, strand):
    """ strand can be either 1 or -1. interval is half open: [forStart, forEnd)
    """
    if strand == 1:
        array[forStart : forEnd] += 1
    else:
        # revStart = chromLength - forEnd
        # revEnd   = chromLength - forStart
        array[totalLength - forEnd : totalLength - forStart] += 1

def pickleData(data, filename, options):
    """ record data to a pickle.
    """
    if filename is None:
        return
    f = open(filename, 'wb')
    cPickle.dump(data, f, 2) # 2 is the format protocol, 2 = binary
    f.close()

def main():
    usage = ('usage: %prog --maf path/to/file.maf --species=species1,species2,...\n\n'
             '%prog is a script that operates on a single maf\n'
             '(multiple alignment format) file and extracts information from\n'
             'the alignments of pairs of sequences. Output is a pickle.\n'
             'The maf sequence name field must have the format\n'
             'species.chromosome\n'
             'e.g.\n'
             'hg19.chr22\n\n'
             'For slightly more detailed examples, check the \n'
             'test.mafCoveragePickleCreator.py unittest.\n\n'
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
    
    alignmentDict = readMaf(options.maf, options)
    pickleData(alignmentDict, options.pickle, options)
    
if __name__ == '__main__':
    main()
