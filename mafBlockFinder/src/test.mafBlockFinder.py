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
import re
import sys
import unittest
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), '../../include/')))
import mafToolsTest as mtt

g_targetSeq = 'target.chr0'
g_header = None
g_headers = ['''##maf version=1 scoring=tba.v8
# tba.v8 (((human chimp) baboon) (mouse rat))

''', 
             '''##maf version=1 scoring=tba.v8
# tba.v8 (((human chimp) baboon) (mouse rat))
''',]

# these are triples of blocks, target positions for True, 
# line numbers where the target occurs.
g_overlappingBlocks = [('''a score=0
# --pos 0
s target.chr0         0 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA
s name.chr1           0 10 +       100 ATGT---ATGCCG
s name2.chr1         50 10 +       100 ATGT---ATGCCG
s name3.chr9         50 10 +       100 ATGTA---TGCCG
s name4.chr&         50 10 +       100 ATG---TATGCCG
s name5 50 10 + 100 ATGTATGCCG

''', 0, [3]),
    ('''a score=0
# --pos 38
s target.chr0        38 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA
s name.chr1           0 10 +       100 ATGT---ATGCCG
s name2.chr1         50 10 +       100 ATGT---ATGCCG
s name3.chr9         50 10 +       100 ATGTA---TGCCG
s name4.chr&         50 10 +       100 ATG---TATGCCG
s name5 50 10 + 100 ATGTATGCCG

''', 38, [3]),
            # the overlap is one base long, right on 50.
            ('''a score=0
# --pos 60
s name                0 10 +       100 ATGTATGC---CG
s name2.chr1         50 10 +       100 ATGTATG---CCG
s name3.chr9         50 10 +       100 ATGTATG---CCG
s name4.chr&         50 10 +       100 ATGTATG---CCG
s name5              50 10 +       100 ATGTATG---CCG
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s target.chr0 158545457 10 - 158545518 ATGTATG---CCG

''', 60, [10]),
            # the overlap is 10 bases long, 51-61
            ('''a score=0
# --pos 70
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

''', 70, [9]),
            # the overlap is 9 bases long, 62-70
            ]
g_nonOverlappingBlocks = [('''a score=23262.0     
s hg18.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon         116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG

''', 0, None),
                          ('''a score=5062.0                    
s hg18.chr7    27699739 6 + 158545518 TAAAGA
s panTro1.chr6 28862317 6 + 161576975 TAAAGA
s baboon         241163 6 +   4622798 TAAAGA 
s mm4.chr6     53303881 6 + 151104725 TAAAGA
s rn3.chr4     81444246 6 + 187371129 taagga

''', 0, None),
                          ('''a score=6636.0
s hg18.chr7    27707221 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA

''', 0, None),
                          ('''a score=0
# --pos 38
s target.chr0        38 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA
s name.chr1           0 10 +       100 ATGT---ATGCCG
s name2.chr1         50 10 +       100 ATGT---ATGCCG
s name3.chr9         50 10 +       100 ATGTA---TGCCG
s name4.chr&         50 10 +       100 ATG---TATGCCG
s name5 50 10 + 100 ATGTATGCCG

''', 37, None),
    ]

def foundLines(lineList, text, g_header):
    f = open(text, 'r')
    n = len(re.findall(re.compile('\n'), g_header))
    for line in f:
        line = line.strip()
        data = line.split(':')
        l = data[0].split(',')[1].split()[1]
        if int(l) - n not in lineList:
            f.close()
            return False
    f.close()
    return True
    
class FindTest(unittest.TestCase):
    def testFind(self):
        """ mafBlockFinder should report information about matching sequences within blocks.
        """
        global g_header
        mtt.makeTempDirParent()
        for i in xrange(0, len(g_overlappingBlocks)):
            tmpDir = os.path.abspath(mtt.makeTempDir('find'))
            testMafPath, g_header = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test.maf')),
                                                 g_overlappingBlocks[i][0], g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = [os.path.abspath(os.path.join(parent, 'test', 'mafBlockFinder'))]
            cmd += ['--maf', testMafPath, '--seq', g_targetSeq, '--pos', '%d' % g_overlappingBlocks[i][1]]
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'found.txt'))]
            mtt.recordCommands([cmd], tmpDir, outPipes=outpipes)
            mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
            self.assertTrue(foundLines(g_overlappingBlocks[i][2], os.path.join(tmpDir, 'found.txt'), g_header))
            mtt.removeDir(tmpDir)
    def testNonFind(self):
        """ mafBlockFinder should not report any lines when blocks do not match.
        """
        global g_header
        mtt.makeTempDirParent()
        for i in xrange(0, len(g_nonOverlappingBlocks)):
            tmpDir = os.path.abspath(mtt.makeTempDir('nonFind'))
            testMafPath, g_header = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test.maf')),
                                                 ''.join(g_nonOverlappingBlocks[i][0]), g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = [os.path.abspath(os.path.join(parent, 'test', 'mafBlockFinder'))]
            cmd += ['--maf', testMafPath, '--seq', g_targetSeq, '--pos', '%d' % g_nonOverlappingBlocks[i][1]]
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'found.txt'))]
            mtt.recordCommands([cmd], tmpDir, outPipes=outpipes)
            mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
            self.assertTrue(mtt.fileIsEmpty(os.path.join(tmpDir, 'found.txt')))
            mtt.removeDir(tmpDir)
    def testMemory1(self):
        """ If valgrind is installed on the system, check for memory related errors (1).
        """
        mtt.makeTempDirParent()
        valgrind = mtt.which('valgrind')
        if valgrind is None:
            return
        global g_header
        for i in xrange(0, len(g_overlappingBlocks)):
            tmpDir = os.path.abspath(mtt.makeTempDir('memory1'))
            testMafPath, g_header = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test.maf')),
                                                 g_overlappingBlocks[i][0], g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = mtt.genericValgrind(tmpDir)
            cmd.append(os.path.abspath(os.path.join(parent, 'test', 'mafBlockFinder')))
            cmd += ['--maf', testMafPath, '--seq', g_targetSeq, '--pos', '%d' % g_overlappingBlocks[i][1]]
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'found.txt'))]
            mtt.recordCommands([cmd], tmpDir, outPipes=outpipes)
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
        global g_header
        for i in xrange(0, len(g_nonOverlappingBlocks)):
            tmpDir = os.path.abspath(mtt.makeTempDir('memory2'))
            testMafPath, g_header = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test.maf')),
                                                 ''.join(g_nonOverlappingBlocks[i][0]), g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = mtt.genericValgrind(tmpDir)
            cmd.append(os.path.abspath(os.path.join(parent, 'test', 'mafBlockFinder')))
            cmd += ['--maf', testMafPath, '--seq', g_targetSeq, '--pos', '%d' % g_nonOverlappingBlocks[i][1]]
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'found.txt'))]
            mtt.recordCommands([cmd], tmpDir, outPipes=outpipes)
            mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
            self.assertTrue(mtt.noMemoryErrors(os.path.join(tmpDir, 'valgrind.xml')))
            mtt.removeDir(tmpDir)

if __name__ == '__main__':
    unittest.main()
