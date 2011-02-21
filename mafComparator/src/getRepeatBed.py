#!/usr/bin/env python
"""

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
