#!/usr/bin/env python

import sys
from sonLib.bioio import fastaRead

if "--help" in sys.argv[1:] or len(sys.argv) == 1:
    print "Script to create a bed file containing the intervals of repeat bases in a  fasta file."
    print "Usage: fastaFile outputBedFile"
    sys.exit(0)

def fn(header):
    return header.split()[0]

def fn2(sequence):
    fn = lambda x : x in [ 'a', 'c', 't', 'g', 'N', 'n']
    i = 0
    while i < len(sequence):
        if fn(sequence[i]):
            j = i+1
            while j<len(sequence) and fn(sequence[j]):
                j+=1
            yield i, j
            i = j
        else:
            i+=1

fileHandle = open(sys.argv[2], 'w')
for sequenceFile in sys.argv[1].split():
    for header, sequence in fastaRead(open(sequenceFile, 'r')):
        sequenceName = fn(header)
        for start, stop, in fn2(sequence):
            fileHandle.write("%s\t%i\t%i\n" % (sequenceName, start, stop))
fileHandle.close()