#!/usr/bin/env python
"""
eval_setDiffDroppedMissing.py
dent earl, dearl (a) soe ucsc edu
3 nov 2010
A script that will take the dropped regions output from mafJoin
(dropped.tab) and the missing pairs from eval_MAFcomparator
(missing.tab) and will look for 'missing pairs' that are not
in the intersection with dropped.tab, i.e. the set difference.
These members outside the intersection represent possible bugs.

"""
##############################
import os
import re
import sys
from optparse import OptionParser

def usage():
    sys.stderr.write('USAGE: %s --droppedFile=fileName --missingFile=fileName \n' % (sys.argv[0]))
    sys.exit(2)

def initOptions(parser):
    parser.add_option('-d', '--droppedFile',dest='droppedFile',
                      help='Regions dropped from maf by mafJoin.')
    parser.add_option('-m', '--missingFile',dest='missingFile',
                      help='Missing pairs from the sampling approach of eval_MAFcomparator.')

def checkOptions(options):
    if (options.droppedFile == None):
        sys.stderr.write('%s: Error, specify --droppedFile.\n' % sys.argv[0])
        usage()
    if (not os.path.exists(options.droppedFile)):
        sys.stderr.write('%s: Error, "%s" does not exist!\n' % (sys.argv[0], options.droppedFile))
        usage()
    options.droppedFile = os.path.abspath(options.droppedFile)

    if (options.missingFile == None):
        sys.stderr.write('%s: Error, specify --missingFile.\n' % sys.argv[0])
        usage()
    if (not os.path.exists(options.missingFile)):
        sys.stderr.write('%s: Error, "%s" does not exist!\n' % (sys.argv[0], options.missingFile))
        usage()
    options.missingFile = os.path.abspath(options.missingFile)

class Drop:
    def __init__(self):
        self.file  = ''
        self.name  = ''
        self.start = -1
        self.end   = -1

def populateDroppedRegions( file ):
    f = open( file )
    dropped = []
    firstLine = True
    for line in f:
        if firstLine:
            firstLine = False
            continue
        line = line.strip()
        t = line.split('\t')
        assert( len(t) == 4 )
        d = Drop()
        d.file  = t[0]
        d.name  = t[1]
        d.start = int( t[2] )
        d.end   = int( t[3] )
        dropped.append(d)
    dropped.sort( key = lambda x: ( x.name, x.start ) )
    return dropped

def printDropped( dropped ):
    for d in dropped:
        print '%s %s %d %d' %( d.file, d.name, d.start, d.end )

class Miss:
    def __init__(self):
        self.file = ''
        self.seq1 = ''
        self.seq2 = ''
        self.pos1 = -1
        self.pos2 = -1

def populateMissingPairs( file ):
    f = open( file )
    pat = re.compile('# Comparing (.*?) to (.*?)$')
    bFile = '' # the second file in the comparison between A and B
    missing = []
    for line in f:
        if line[0] == '#':
            r = re.match(pat, line)
            if r:
                bFile = r.group(2)
            continue
        line = line.strip()
        t = line.split('\t')
        m = Miss()
        m.file = bFile
        m.seq1 = t[0]
        m.pos1 = int( t[1] )
        m.seq2 = t[2]
        m.pos2 = int( t[3] )
        missing.append( m )
    missing.sort( key = lambda x: ( x.file, x.seq1, x.pos1 ) )
    return missing

def printMissing( missing ):
    for m in missing:
        print '%s %d %s %d' %( m.seq1, m.pos1, m.seq2, m.pos2 )

def doSetDifference( dropped, missing ):
    inter = intersect( dropped, missing )
    intersectC = setDiff( dropped, missing, inter )
    return intersectC

def setDiff( dropped, missing, intersect ):
    intersectDict = {}
    diff = []
    for i in intersect:
        intersectDict[ i ] = True
    for m in missing:
        if m.seq1 == m.seq2:
            continue
        if m not in intersect:
            diff.append( m )
    return diff

# hey, check out bisect
# wait, create a Dropped dict keyed with names and valued with an array
# 

def intersect( dropped, missing ):
    intersect = []
    for m in missing:
        if m.seq1 == m.seq2:
            continue
        for d in dropped:
            if m.seq1 != d.name and m.seq2!= d.name:
                continue
            if m.seq1 == d.name:
                if m.pos1 < d.start:
                    continue
                if m.pos1 >= d.end:
                    continue
                intersect.append( m )
            if m.seq2 == d.name:
                if m.pos2 < d.start:
                    continue
                if m.pos2 >= d.end:
                    continue
                intersect.append( m )
    return intersect
    
def main():
    parser=OptionParser()
    initOptions(parser)
    (options, args) = parser.parse_args()
    checkOptions(options)

    dropped = []
    dropped = populateDroppedRegions( options.droppedFile )
    #printDropped( dropped )

    missing = []
    missing = populateMissingPairs( options.missingFile )
    #printMissing( missing )

    setDifference = []
    setDifference = doSetDifference( dropped, missing )
    printMissing( setDifference )

if __name__ == "__main__":
    main()
