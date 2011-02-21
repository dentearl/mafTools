#!/usr/bin/env python
"""
setDiffDroppedMissing.py
dent earl, dearl (a) soe ucsc edu
3 nov 2010
A script that will take the dropped regions output from mafJoin
(dropped.tab) and the missing pairs from MAFcomparator
(missing.tab) and will look for 'missing pairs' that are not
in the intersection with dropped.tab, i.e. the set difference.
These members outside the intersection represent possible bugs.


Copyright (C) 2009-2011 by 
Dent Earl (dearl@soe.ucsc.edu, dentearl@gmail.com)
Benedict Paten (benedict@soe.ucsc.edu, benedictpaten@gmail.com)
Mark Diekhans (markd@soe.ucsc.edu)
... and other members of the Reconstruction Team of David Haussler's 
lab (BME Dept. UCSC).

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.


"""
##############################
import bisect
import os
import re
import sys
from optparse import OptionParser

def usage():
    sys.stderr.write('USAGE: %s --droppedFile=fileName --missingFile=fileName \n' % (sys.argv[0]))
    sys.exit(2)

def initOptions(parser):
    parser.add_option('-d', '--droppedFile', dest='droppedFile',
                      help='Regions dropped from maf by mafJoin.')
    parser.add_option('-m', '--missingFile', dest='missingFile',
                      help='Missing pairs from the sampling approach of MAFcomparator.')
    parser.add_option('-v', '--verbose', dest='isVerbose', action='store_true',
                      default=False, help='Prints out extra information to STDERR.')

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
#    def __eq__( self, other ):
#        return self.start == other
    def __lt__( self, other ):
        return self.start < other
    def __gt__( self, other ):
        return self.start > other
    def __le__( self, other ):
        return self.__eq__( other ) or self.__lt__( other )
    def __ge__( self, other ):
        return self.__eq__( other ) or self.__gt__( other )

def populateDroppedRegions( file, isVerbose ):
    """ populateDroppedRegions produces a dict keyed by sequence name and valued
    with an array that contains Drop() objects.
    """
    f = open( file )
    dropped = {}
    for line in f:
        t = line.split('\t')
        assert( len(t) == 4 )
        if t[0] == 'maf':
            continue
        d = Drop()
        d.file  = t[0]
        d.name  = t[1]
        d.start = int( t[2] )
        d.end   = int( t[3] )
        if not d.name in dropped:
            dropped[ d.name ] = []
        dropped[ d.name ].append( d )
    for n in dropped:
        dropped[ n ].sort()
    return dropped

def printDropped( dropped ):
    print '#Dropped from\tname\tstart\tend'
    for n in dropped:
        for d in dropped[ n ]:
            print '%s\t%s\t%d\t%d' %( d.file, d.name, d.start, d.end )

class Miss:
    def __init__(self):
        self.file = ''
        self.seq1 = ''
        self.seq2 = ''
        self.pos1 = -1
        self.pos2 = -1

def populateMissingPairs( file, isVerbose ):
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
        t = line.split('\t')
        m = Miss()
        m.file = bFile
        m.seq1 = t[0]
        m.pos1 = int( t[1] )
        m.seq2 = t[2]
        m.pos2 = int( t[3] )
        if m.seq1 != m.seq2: # for the time being we ignore self-self alignments
            missing.append( m )
    missing.sort( key= lambda x: ( x.file, x.seq1, x.pos1 ) )
    return missing

def printMissing( missing ):
    print '#Missing from\tseq1\tpos1\tseq2\tpos2'
    for m in missing:
        print '%s\t%s\t%d\t%s\t%d' %( m.file, m.seq1, m.pos1, m.seq2, m.pos2 )

def doSetDiff( dropped, missing, isVerbose ):
    message(isVerbose, 'starting intersection... ')
    intersectDict = intersect( dropped, missing )
    message(isVerbose, ' OK\n')
    message(isVerbose, 'starting set difference... ')
    setDifference = setDiff( missing, intersectDict )
    message(isVerbose, ' OK\n')
    return setDifference

def setDiff( missing, intersect ):
    diff = []
    for m in missing:
        if m not in intersect:
            diff.append( m )
    return diff

def find_le(a, x):
    'Find rightmost value less than or equal to x'
    i = bisect.bisect_right(a, x)
    if i:
        return i-1
    return 0

def intersect( dropped, missing ):
    """ intersect() takes the hash-list structure of dropped regions
    and the list of missing pairs and then uses the bisect module
    to search for missing pairs in the dropped regions.
    """
    intersect = {}
    for m in missing:
        if m.seq1 in dropped:
            i = find_le( dropped[ m.seq1 ], m.pos1 )
            if m.pos1 < dropped[ m.seq1 ][ i ].start:
                continue
            if m.pos1 < dropped[ m.seq1 ][ i ].end:
                intersect[ m ] = True
        if m.seq2 in dropped:
            i = find_le( dropped[ m.seq2 ], m.pos2 )
            if m.pos2 < dropped[ m.seq2 ][ i ].start:
                continue
            if m.pos2 < dropped[ m.seq2 ][ i ].end:
                intersect[ m ] = True
    return intersect

def message(isVerbose, s):
    if isVerbose:
        sys.stderr.write(s)
        sys.stderr.flush()
    
def main():
    parser=OptionParser()
    initOptions(parser)
    (options, args) = parser.parse_args()
    checkOptions(options)
    message(options.isVerbose, 'loading \'dropped\' file... ')
    dropped = {}
    dropped = populateDroppedRegions( options.droppedFile, options.isVerbose )
    message(options.isVerbose, 'OK\n')
    # printDropped( dropped )

    message(options.isVerbose, 'loading \'missing\' file... ')
    missing = []
    missing = populateMissingPairs( options.missingFile, options.isVerbose )
    message(options.isVerbose, 'OK\n')
    # printMissing( missing )

    setDifference = []
    setDifference = doSetDiff( dropped, missing, options.isVerbose )
    printMissing( setDifference )

if __name__ == "__main__":
    main()
