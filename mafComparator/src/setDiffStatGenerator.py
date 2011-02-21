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
##############################
import os
import re
import sys
from optparse import OptionParser

def initOptions(parser):
    parser.add_option('-s', '--second', dest='second', action='store_true',
                      default=False, help='Insead of taking first sequence in the pair, take the second.')

def checkOptions(options):
    pass

def readStream(options):
    sizeDict = {}
    positionDict = {}
    nameDict = {}
    offset = 0
    if options.second:
        offset = 2
    startPos = -1
    currLength  = 0
    previousPos = -1
    for line in sys.stdin:
        t = line.split('\t')
        if t[ 0 ][0] == '#':
            continue
        if int( t[ 2 + offset ] ) == previousPos:
            continue
        if int( t[ 2 + offset ] ) == previousPos + 1:
            currLength += 1
            previousPos = int( t[ 2 + offset ] )
            continue
        if currLength in sizeDict:
            sizeDict[ currLength ] += 1
            positionDict[ currLength ] = startPos
            nameDict[ currLength ] = t[ 1 + offset ]
        else:
            sizeDict[ currLength ]  = 1
            positionDict[ currLength ] = startPos
            nameDict[ currLength ] = t[ 1 + offset ]
        currLength = 1
        startPos = int( t[ 2 + offset ] )
        previousPos = int( t[ 2 + offset ] )
    if currLength in sizeDict:
        sizeDict[ currLength ] += 1
    else:
        sizeDict[ currLength ]  = 1
    keys = sizeDict.keys()
    keys.sort()
    print ' bp\tCount\tMR pos\tMR name'
    for k in keys:
        print '%3d\t%d\t%d\t%s' % ( k, sizeDict[ k ], positionDict[ k ], nameDict[ k ] )

def main():
    parser=OptionParser()
    initOptions(parser)
    (options, args) = parser.parse_args()
    checkOptions(options)

    readStream( options )


if __name__ == "__main__":
    main()
