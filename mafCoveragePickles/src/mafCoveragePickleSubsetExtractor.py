#!/usr/bin/env python2.6
"""
mafCoveragePickleSubsetExtractor.py
10 January 2012
dent earl

These pickles are large and in order to plot from one
you need to read the entire thing. So, this script allows
you to read one of the large pickles and pull out one array
and store it in its own pickle to make plotting faster.

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
import mafIndelDistribution as mid
import numpy
from optparse import OptionParser
import os

def initOptions(parser):
    parser.add_option('--pickle', dest = 'pickle',
                      help = 'location of pickle to read')
    parser.add_option('--out', dest = 'out',
                      help = 'location of pickle to write')
    parser.add_option('--spp1', dest = 'spp1',
                      help = 'Species 1.')
    parser.add_option('--spp1chr', dest = 'spp1chr',
                      help = 'Species 1 chromosome.')
    parser.add_option('--spp2', dest = 'spp2',
                     help = 'Species 2.')

def checkOptions(options, args, parser):
    for k, v in [('spp1', options.spp1), ('spp2', options.spp2),
                 ('spp1chr', options.spp1chr), ('pickle', options.pickle),
                 ('out', options.out)]:
        if v is None:
            parser.error('specify --%s' % k)
    if not os.path.exists(options.pickle):
        parser.error('--pickle %s does not exist.' % options.pickle)

def main():
    usage = ('usage: %prog --pickle path/to/maf.pickle --spp1 species1 --spp1chr chrom1 '
             '--spp2 species2 --out subset.pickle \n\n'
             '%prog opens a pickled maf and pulls out one species+chromosome to other species\n'
             'alignment array and stores it in a new pickle. The total pickles can get quite\n'
             'big and having to read through the entire thing when operating or plotting only\n'
             'a small part of pickle is annoying. Instead use this script to pull out the part\n'
             'you\'re working on and then iterate on that subset.')
    parser = OptionParser(usage = usage)
    initOptions(parser)
    options, args = parser.parse_args()
    checkOptions(options, args, parser)
    
    dataDict = {options.pickle : mcpp.readPickle(options.pickle)}
    speciesSet, speciesChrDict = mcpp.parseDataDict(dataDict, options, parser)
    
    export = {options.spp1 : 
              {options.spp1chr :
                   {options.spp2 : dataDict[options.pickle][options.spp1][options.spp1chr][options.spp2]}}}
    mid.pickleData(export, options.out, options)
    
if __name__ == '__main__':
    main()
