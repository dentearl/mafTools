#!/usr/bin/env python2.7
"""
mafCoveragePicklePlotter.py
6 January 2012
dent earl

script to create proportional coverage plots of one 
species onto another species+chromosome using
alignment pickles.

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
import matplotlib
matplotlib.use('Agg')
import cPickle
from math import ceil
from mafCoveragePickleCreator import readPickles
import matplotlib.backends.backend_pdf as pltBack
import matplotlib.pyplot as plt
import numpy
from optparse import OptionParser
import os

def initOptions(parser):
    parser.add_option('--spp1', dest='spp1',
                      help='Species 1.')
    parser.add_option('--spp1chr', dest='spp1chr',
                      help='Species 1 chromosome.')
    parser.add_option('--spp2', dest='spp2',
                     help='Species 2.')
    parser.add_option('--windowLength', dest='windowLength', 
                      default=500 * 10**3, type='int',
                      help='Size of sliding window. default=%default')
    parser.add_option('--windowOverlap', dest='windowOverlap', 
                      default=250 * 10**3, type='int',
                      help='Size of sliding window step. default=%default')
    parser.add_option('--dpi', dest='dpi', default=300,
                      type='int',
                      help='Dots per inch of the output, if --outFormat is all or png. default=%default')
    parser.add_option('--outFormat', dest='outFormat', default='pdf',
                      type='string',
                      help='output format [pdf|png|all|eps]. default=%default')
    parser.add_option('--out', dest='out', default='myPlot',
                      type='string',
                      help='filename where figure will be created. No extension needed. default=%default')
    
def checkOptions(options, args, parser):
    for k, v in [('spp1', options.spp1), ('spp2', options.spp2),
                 ('spp1chr', options.spp1chr)]:
        if v is None:
            parser.error('specify --%s' % k)
    if options.dpi < 72:
        parser.error('--dpi %d too low. Must be at least 72.' % options.dpi)
    if options.outFormat not in ('pdf', 'png', 'eps', 'all'):
        parser.error('Unrecognized output format: %s. Choose one from: pdf png eps all.' % options.outFormat)
    if (options.out.endswith('.png') or options.out.endswith('.pdf') or 
         options.out.endswith('.eps')):
        options.out = options.out[:-4]
    options.windowStep = options.windowLength - options.windowOverlap
    for a in args:
        if not os.path.exists(a):
            parser.error('File %s does not exist.' % a)

def parseDataDict(dataDict, options, parser):
    """ look in the data and find out what species are available
    and what chromosomes they possess. This will help with input 
    options checking.
    """
    speciesSet = set()
    speciesChrDict = {}
    
    for p in dataDict:
        for s1 in dataDict[p]:
            if s1 not in speciesSet:
                speciesSet.add(s1)
            if s1 not in speciesChrDict:
                speciesChrDict[s1] = set()
            for c in dataDict[p][s1]:
                if c not in speciesChrDict[s1]:
                    speciesChrDict[s1].add(c)
                for s2 in dataDict[p][s1][c]:
                    if s2 not in speciesSet:
                        speciesSet.add(s2)
    for u,v in [('spp1', options.spp1), ('spp2', options.spp2)]:
        if v not in speciesSet:
            parser.error('--%s %s is not found any of the input pickles' % (u, v))
    if options.spp1chr not in speciesChrDict[options.spp1]:
        parser.error('--spp1chr %s is not associated with species %s in any input pickles' 
                     % (options.spp1chr, options.spp1))
    return speciesSet, speciesChrDict

def initImage(width, height, options):
    """
    initImage takes a width and height and returns
    both a fig and pdf object. options must contain outFormat,
    and dpi
    """
    pdf = None
    if options.outFormat == 'pdf' or options.outFormat == 'all':
        pdf = pltBack.PdfPages(options.out + '.pdf')
    fig = plt.figure(figsize=(width, height), dpi=options.dpi, facecolor='w')
    return (fig, pdf)

def writeImage(fig, pdf, options):
    """
    writeImage assumes options contains outFormat and dpi.
    """
    if options.outFormat == 'pdf':
        fig.savefig(pdf, format='pdf')
        pdf.close()
    elif options.outFormat == 'png':
        fig.savefig(options.out + '.png', format='png', dpi=options.dpi)
    elif options.outFormat == 'all':
        fig.savefig(pdf, format='pdf')
        pdf.close()
        fig.savefig(options.out + '.png', format='png', dpi=options.dpi)
        fig.savefig(options.out + '.eps', format='eps')
    elif options.outFormat == 'eps':
        fig.savefig(options.out + '.eps', format='eps')

def establishAxis(fig, options):
   """ creates the axis that the plot is drawn upon
   """
   options.axLeft  = 0.09
   options.axWidth = 0.88
   options.axBottom = 0.12
   options.axHeight = 0.80
   ax = fig.add_axes([options.axLeft, options.axBottom,
                       options.axWidth, options.axHeight])
   return ax

def drawAnnotations(ax, plotList, filenames, options):
    ax.set_title('%s onto %s %s' % (options.spp2, options.spp1, options.spp1chr))
    plt.ylabel('Prop. Coverage (Win. = %s Overl. = %s)' % (prettyBase(options.windowLength, 1), 
                                                           prettyBase(options.windowOverlap, 1)))
    plt.xlabel('Position along %s %s (Mb)' % (options.spp1, options.spp1chr))

    if len(plotList) > 0:
        leg = plt.legend(plotList, filenames)
        plt.setp(leg.get_texts(), fontsize='small') # legend fontsize
        leg._drawFrame = False

def drawOne(xdata, data, color, plotList, ax, options):
    xdata /= 10**6
    plotList.append(ax.plot(xdata, data, color=color))

def prettyBase(n, decimals = 2):
    """ takes an int, assumed to be a number of basepairs and 
    returns a string with proper units
    """
    if n >= 10**9:
        return '%.*f Gb' % (decimals, n / 10.0**9)
    elif n >= 10**6:
        return '%.*f Mb' % (decimals, n / 10.0**6)
    elif n >= 10**3:
        return '%.*f Kb' % (decimals, n / 10.0**3)
    else:
        return '%d' % n

def window(a, options):
    """ take a numpy array, `a', and return a numpy array that
    is the result of a sliding window pass on a.
    Note that this method of sliding window will result in edge effects
    on the right end of the array. These will not be noticeable as long as the
    window size / total length(a) is small, say less than 0.02
    """
    aLength = len(a)
    wind = numpy.zeros(ceil(aLength/float(options.windowStep)), dtype=numpy.float32)
    xdata = numpy.zeros(ceil(aLength/float(options.windowStep)), dtype=numpy.float32)
    halfWind = (options.windowLength / 2.0)
    i = -1
    for start in xrange(0, aLength, options.windowStep):
        i += 1
        stop = min(start + options.windowStep, aLength - 1)
        wind[i] = sum(a[start:stop]) / float(stop - start)
        xdata[i] = halfWind + i * options.windowStep
    return xdata, wind

def main():
    usage = ('usage: %prog --spp1 species1 --spp1chr chrom1 --spp2 species2 maf1.pickle maf2.pickle ...\n\n'
             '%prog takes a species+chromosome to other species alignment pair and a set of \n'
             'alignment pickles (as positional arguments) and produces a plot showing the\n'
             'proportional coverage of --spp2 onto --spp1 along --spp1chr.')
    parser = OptionParser(usage=usage)
    initOptions(parser)
    options, args = parser.parse_args()
    checkOptions(options, args, parser)
    
    dataDict = readPickles(args)
    speciesSet, speciesChrDict = parseDataDict(dataDict, options, parser)
    
    fig, pdf = initImage(8.5, 5.0, options)
    ax = establishAxis(fig, options)
    plotList = []
    filenames = []
    # colorList = ['#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78']
    colorList = ['#4B4C5E', # dark slate gray
                 '#9edae5', # light blue 
                 '#7F80AB', # purple-ish slate blue
                 '#c7c7c7', # light gray
                 '#ff7f0e', # bright orange
                 '#ffbb78', # light orange
                 '#9467bd', # dark purple
                 '#c5b0d5',  # light purple
                 ]
    for i, a in enumerate(args, 0):
        if options.spp2 not in dataDict[a][options.spp1][options.spp1chr]:
            print 'Warning: %s onto %s %s is not present in %s' % (options.spp2, options.spp1,
                                                                   options.spp1chr, a)
            continue
        filenames.append(a)
        xdata, ydata = window(dataDict[a][options.spp1][options.spp1chr][options.spp2], options)
        drawOne(xdata, ydata, colorList[i], plotList, ax, options)
    drawAnnotations(ax, plotList, filenames, options)
    writeImage(fig, pdf, options)

if __name__ == '__main__':
    main()
