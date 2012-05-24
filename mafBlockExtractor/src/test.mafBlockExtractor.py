##################################################
# Copyright (C) 2012 by 
# Dent Earl (dearl@soe.ucsc.edu, dentearl@gmail.com)
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
##################################################
import os
import random
import sys
import unittest
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), '../../include/')))
import mafToolsTest as mtt

g_targetSeq = 'target.chr0'
g_targetRange = (50, 70) # zero based, inclusive

g_headers = ['''##maf version=1 scoring=tba.v8
# tba.v8 (((human chimp) baboon) (mouse rat))

''',
             '''##maf version=1 scoring=tba.v8
# tba.v8 (((human chimp) baboon) (mouse rat))
''',
             '''##maf version=1

''',]

g_overlappingBlocks = ['''a score=0
s target.chr0        38 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA
s name.chr1           0 10 +       100 ATGT---ATGCCG
s name2.chr1         50 10 +       100 ATGT---ATGCCG
s name3.chr9         50 10 +       100 ATGTA---TGCCG
s name4.chr&         50 10 +       100 ATG---TATGCCG
s name5 50 10 + 100 ATGTATGCCG

''',
            # the overlap is one base long, right on 50.
            '''a score=0
s name                0 10 +       100 ATGTATGC---CG
s name2.chr1         50 10 +       100 ATGTATG---CCG
s name3.chr9         50 10 +       100 ATGTATG---CCG
s name4.chr&         50 10 +       100 ATGTATG---CCG
s name5              50 10 +       100 ATGTATG---CCG
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s target.chr0 158545457 10 - 158545518 ATGTATG---CCG

''',
            # the overlap is 10 bases long, 51-61
            '''a score=0
s name               10 10 +       100 ATGTAT---GCCG
s name2.chr1         50 10 +       100 ATGTAT---GCCG
s name3.chr9         50 10 +       100 ATGTAT---GCCG
s name4.chr&         50 10 +       100 ATGTAT---GCCG
s name5              50 10 +       100 ATGTAT---GCCG
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s target.chr0        62  9 + 158545518 gca---gaa-aca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s hg16.chr7    27707221 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA

''',
            # the overlap is 9 bases long, 62-70
            ]
g_nonOverlappingBlocks = ['''a score=23262.0     
s hg18.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon         116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG

''',
                          '''a score=5062.0                    
s hg18.chr7    27699739 6 + 158545518 TAAAGA
s panTro1.chr6 28862317 6 + 161576975 TAAAGA
s baboon         241163 6 +   4622798 TAAAGA 
s mm4.chr6     53303881 6 + 151104725 TAAAGA
s rn3.chr4     81444246 6 + 187371129 taagga

''',
                          '''a score=6636.0
s hg18.chr7    27707221 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA

''',
    ]

def mafIsExtracted(maf):
    f = open(maf)
    lastLine = mtt.processHeader(f)
    for i in xrange(0, len(g_overlappingBlocks)):
        b = mtt.extractBlockStr(f, lastLine)
        lastLine = None
        if b not in g_overlappingBlocks:
            print 'dang'
            print b
            print '!='
            print g_overlappingBlocks
            return False
    return True
    
class ExtractionTest(unittest.TestCase):
    def testExtraction(self):
        """ mafBlockExtractor should output blocks that meet the criteria for extraction. That is they contain the target sequence and have at least one base in the target range.
        """
        mtt.makeTempDirParent()
        for i in xrange(0, 10):
            shuffledBlocks = []
            tmpDir = os.path.abspath(mtt.makeTempDir('extraction'))
            order = [1] * len(g_overlappingBlocks) + [0] * len(g_nonOverlappingBlocks)
            random.shuffle(order)
            random.shuffle(g_overlappingBlocks)
            random.shuffle(g_nonOverlappingBlocks)
            j, k = 0, 0
            for isOverlapping in order:
                if isOverlapping:
                    shuffledBlocks.append(g_overlappingBlocks[j])
                    j += 1
                else:
                    shuffledBlocks.append(g_nonOverlappingBlocks[k])
                    k += 1
            testMaf = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test.maf')),
                                   ''.join(shuffledBlocks), g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = [os.path.abspath(os.path.join(parent, 'test', 'mafBlockExtractor'))]
            cmd += ['--maf', os.path.abspath(os.path.join(tmpDir, 'test.maf')),
                    '--seq', g_targetSeq, '--start', '%d' % g_targetRange[0], 
                    '--stop', '%d' % g_targetRange[1]]
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'extracted.maf'))]
            mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
            self.assertTrue(mafIsExtracted(os.path.join(tmpDir, 'extracted.maf')))
            mtt.removeDir(tmpDir)
    def testNonExtraction(self):
        """ mafBlockExtractor should not extract blocks when they do not match.
        """
        mtt.makeTempDirParent()
        for i in xrange(0, 10):
            tmpDir = os.path.abspath(mtt.makeTempDir())
            random.shuffle(g_nonOverlappingBlocks)
            testMaf = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test.maf')),
                                   ''.join(g_nonOverlappingBlocks), g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = [os.path.abspath(os.path.join(parent, 'test', 'mafBlockExtractor'))]
            cmd += ['--maf', os.path.abspath(os.path.join(tmpDir, 'test.maf')),
                    '--seq', g_targetSeq, '--start', '%d' % g_targetRange[0], 
                    '--stop', '%d' % g_targetRange[1]]
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'extracted.maf'))]
            mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
            self.assertTrue(mtt.fileIsEmpty(os.path.join(tmpDir, 'extracted.maf')))
            mtt.removeDir(tmpDir)
    def testMemory1(self):
        """ If valgrind is installed on the system, check for memory related errors (1).
        """
        mtt.makeTempDirParent()
        valgrind = mtt.which('valgrind')
        if valgrind is None:
            return
        for i in xrange(0, 10):
            shuffledBlocks = []
            tmpDir = os.path.abspath(mtt.makeTempDir())
            order = [1] * len(g_overlappingBlocks) + [0] * len(g_nonOverlappingBlocks)
            random.shuffle(order)
            random.shuffle(g_overlappingBlocks)
            random.shuffle(g_nonOverlappingBlocks)
            j, k = 0, 0
            for isOverlapping in order:
                if isOverlapping:
                    shuffledBlocks.append(g_overlappingBlocks[j])
                    j += 1
                else:
                    shuffledBlocks.append(g_nonOverlappingBlocks[k])
                    k += 1
            testMaf = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test.maf')),
                                   ''.join(shuffledBlocks), g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = mtt.genericValgrind(tmpDir)
            cmd.append(os.path.abspath(os.path.join(parent, 'test', 'mafBlockExtractor')))
            cmd += ['--maf', os.path.abspath(os.path.join(tmpDir, 'test.maf')),
                    '--seq', g_targetSeq, '--start', '%d' % g_targetRange[0], 
                    '--stop', '%d' % g_targetRange[1]]
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'extracted.maf'))]
            mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
            self.assertTrue(mtt.noMemoryErrors(os.path.join(tmpDir, 'valgrind.xml')))
            mtt.removeDir(tmpDir)
    def testMemory2(self):
        """ If valgrind is installed on the system, check for memory related errors (2).
        """
        mtt.makeTempDirParent()
        valgrind = mtt.which('valgrind')
        if valgrind is None:
            return
        for i in xrange(0, 10):
            tmpDir = os.path.abspath(mtt.makeTempDir())
            random.shuffle(g_nonOverlappingBlocks)
            testMaf = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test.maf')),
                                   ''.join(g_nonOverlappingBlocks), g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = mtt.genericValgrind(tmpDir)
            cmd.append(os.path.abspath(os.path.join(parent, 'test', 'mafBlockExtractor')))
            cmd += ['--maf', os.path.abspath(os.path.join(tmpDir, 'test.maf')),
                    '--seq', g_targetSeq, '--start', '%d' % g_targetRange[0], 
                    '--stop', '%d' % g_targetRange[1]]
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'extracted.maf'))]
            mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
            self.assertTrue(mtt.noMemoryErrors(os.path.join(tmpDir, 'valgrind.xml')))
            mtt.removeDir(tmpDir)

if __name__ == '__main__':
    unittest.main()
