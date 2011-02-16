#!/usr/bin/env python
"""

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
