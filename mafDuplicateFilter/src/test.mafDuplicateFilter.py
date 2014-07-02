##################################################
# Copyright (C) 2013 by 
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
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), '../../lib/')))
import mafToolsTest as mtt

g_headers = ['''##maf version=1 scoring=tba.v8
# tba.v8 (((human chimp) baboon) (mouse rat))

''',
             '''##maf version=1 scoring=tba.v8
# tba.v8 (((human chimp) baboon) (mouse rat))
''']

g_duplicateBlocks = [('''a score=0
#dup block 1, name4 is duplicate
s target.chr0        38 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 gcagctgaaaaca
s name.chr1           0 13 +       100 gcagctgaaaaca
s name2.chr1         50 13 +       100 gcagctgaaaaca
s name3.chr9         50 13 +       100 gcagctgaaaaca
s name4.chr&         50 13 +       100 gcagctgaaaaca
s name4.chrA         50 13 +       100 gcagctgaaaacT

''',
                      '''a score=0
#dup block 1, name4 is duplicate
s target.chr0        38 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 gcagctgaaaaca
s name.chr1           0 13 +       100 gcagctgaaaaca
s name2.chr1         50 13 +       100 gcagctgaaaaca
s name3.chr9         50 13 +       100 gcagctgaaaaca
s name4.chr&         50 13 +       100 gcagctgaaaaca

''',),
                     ('''a score=0
#dup block 2, target is duplicate
s name                0 13 +       100 gcagctgaaaaca
s name2.chr1         50 13 +       100 gcagctgaaaaca
s name3.chr9         50 13 +       100 gcagctgaaaaca
s name4.chr&         50 13 +       100 gcagctgaaaaca
s name5              50 13 +       100 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s target.chr0 158545457 13 - 158545518 gcagctgaaaacT
s target.chr1 158545457 13 - 158545518 gcagctgaaaaca

''',
                      '''a score=0
#dup block 2, target is duplicate
s name                0 13 +       100 gcagctgaaaaca
s name2.chr1         50 13 +       100 gcagctgaaaaca
s name3.chr9         50 13 +       100 gcagctgaaaaca
s name4.chr&         50 13 +       100 gcagctgaaaaca
s name5              50 13 +       100 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s target.chr1 158545457 13 - 158545518 gcagctgaaaaca

''',),
                     ('''a score=0
#dup block 3, panTro1 and baboon are duplicates
s name               10 13 +       100 gcagctgaaaaca
s name2.chr1         50 13 +       100 gcagctgaaaaca
s name3.chr9         50 13 +       100 gcagctgaaaaca
s name4.chr&         50 13 +       100 gcagctgaaaaca
s name5              50 13 +       100 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 Acagctgaatact
s target.chr0        62  9 + 158545518 gca---gaa-aca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s hg16.chr7    27707221 13 + 158545518 gcagctgaaaaca
s panTro1.chr7 28869787 13 + 161576975 gcagctgaatact
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA

''',
                      '''a score=0
#dup block 3, panTro1 and baboon are duplicates
s name               10 13 +       100 gcagctgaaaaca
s name2.chr1         50 13 +       100 gcagctgaaaaca
s name3.chr9         50 13 +       100 gcagctgaaaaca
s name4.chr&         50 13 +       100 gcagctgaaaaca
s name5              50 13 +       100 gcagctgaaaaca
s target.chr0        62  9 + 158545518 gca---gaa-aca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s hg16.chr7    27707221 13 + 158545518 gcagctgaaaaca
s panTro1.chr7 28869787 13 + 161576975 gcagctgaatact
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA

''',),
                     ('''a score=0
#dup block 4, name, panTro1 and baboon are duplicates
s name               10 13 +       100 gcagctgaaaaca
s name.chr1          50 13 +       100 gcagctgaaaact
s name.chr2          50 13 +       100 gcagctgaaaact
s name3.chr9         50 13 +       100 gcagctgaaaaca
s name4.chr&         50 13 +       100 gcagctgaaaaca
s name5              50 13 +       100 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 Acagctgaatact
s target.chr0        62  9 + 158545518 gca---gaa-aca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s hg16.chr7    27707221 13 + 158545518 gcagctgaaaaca
s panTro1.chr7 28869787 13 + 161576975 gcagctgaatacT
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA

''',
                      '''a score=0
#dup block 4, name, panTro1 and baboon are duplicates
s name               10 13 +       100 gcagctgaaaaca
s name3.chr9         50 13 +       100 gcagctgaaaaca
s name4.chr&         50 13 +       100 gcagctgaaaaca
s name5              50 13 +       100 gcagctgaaaaca
s target.chr0        62  9 + 158545518 gca---gaa-aca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s hg16.chr7    27707221 13 + 158545518 gcagctgaaaaca
s panTro1.chr7 28869787 13 + 161576975 gcagctgaatacT
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA

''',),
                     ('''a score=0
#dup block 1, name4 is duplicate
s target.chr0        38 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 gcagctgannnnn
s name.chr1           0 13 +       100 gcagctgaaaacN
s name2.chr1         50 13 +       100 gcagctgaaaacN
s name3.chr9         50 13 +       100 gcagctgaaaaca
s name4.chr&         50 13 +       100 gcagctgaaaaca
s name4.chrA         50 13 +       100 gcagctgaaaacT

''',
                      '''a score=0
#dup block 1, name4 is duplicate
s target.chr0        38 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 gcagctgannnnn
s name.chr1           0 13 +       100 gcagctgaaaacN
s name2.chr1         50 13 +       100 gcagctgaaaacN
s name3.chr9         50 13 +       100 gcagctgaaaaca
s name4.chr&         50 13 +       100 gcagctgaaaaca

''',),
            ]
g_nonDuplicateBlocks = ['''a score=23262.0
#non-dup block 1
s hg18.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon         116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG

''',
                          '''a score=5062.0
#non-dup block 2
s hg18.chr7    27699739 6 + 158545518 TAAAGA
s panTro1.chr6 28862317 6 + 161576975 TAAAGA
s baboon         241163 6 +   4622798 TAAAGA
# ignore this comment line
s mm4.chr6     53303881 6 + 151104725 TAAAGA
s rn3.chr4     81444246 6 + 187371129 taagga
q i dont remember what q lines do, but we should be ignoring them.

''',
                          '''a score=6636.0
#non-dup block 3
# this comment line should not screw anything up
s hg18.chr7    27707221 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA
# nor should this comment line.

''',
    ]

def mafIsFiltered(maf, blockList):
    f = open(maf)
    lastLine = mtt.processHeader(f)
    for i in xrange(0, len(blockList)):
        # walk through the maf, assessing the equivalence to the blockList items
        b = mtt.extractBlockStr(f, lastLine)
        lastLine = None
        if b != blockList[i]:
            print 'dang'
            print 'observed:'
            print b
            print '!='
            print 'expected:'
            print blockList[i]
            return False
    return True
class DuplicationFilterTest(unittest.TestCase):
    def testFilter(self):
        """ mafDuplicateFilter should filter out duplicates in blocks according to sequence similarity to the consensus.
        """
        mtt.makeTempDirParent()
        for i in xrange(0, 10):
            shuffledBlocks = []
            expectedOutput = []
            tmpDir = os.path.abspath(mtt.makeTempDir('filter'))
            order = [1] * len(g_duplicateBlocks) + [0] * len(g_nonDuplicateBlocks)
            random.shuffle(order)
            random.shuffle(g_duplicateBlocks)
            random.shuffle(g_nonDuplicateBlocks)
            j, k = 0, 0
            for dupBlock in order:
                if dupBlock:
                    shuffledBlocks.append(g_duplicateBlocks[j][0])
                    expectedOutput.append(g_duplicateBlocks[j][1])
                    j += 1
                else:
                    shuffledBlocks.append(g_nonDuplicateBlocks[k])
                    expectedOutput.append(g_nonDuplicateBlocks[k])
                    k += 1
            testMaf = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test.maf')),
                                   ''.join(shuffledBlocks), g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = [os.path.abspath(os.path.join(parent, 'test', 'mafDuplicateFilter')), 
                   '--maf', os.path.abspath(os.path.join(tmpDir, 'test.maf'))]
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'filtered.maf'))]
            mtt.recordCommands([cmd], tmpDir, outPipes=outpipes)
            mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
            self.assertTrue(mafIsFiltered(os.path.join(tmpDir, 'filtered.maf'), expectedOutput))
            mtt.removeDir(tmpDir)
    def testNonFilter(self):
        """ mafDuplicateFilter should not filter out any sequences from blocks when there are no duplicates.
        """
        mtt.makeTempDirParent()
        for i in xrange(0, 10):
            tmpDir = os.path.abspath(mtt.makeTempDir('nonFilter'))
            random.shuffle(g_nonDuplicateBlocks)
            testMaf = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test.maf')),
                                   ''.join(g_nonDuplicateBlocks), g_headers)
            expectedOutput = g_nonDuplicateBlocks
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = [os.path.abspath(os.path.join(parent, 'test', 'mafDuplicateFilter')), 
                   '--maf', os.path.abspath(os.path.join(tmpDir, 'test.maf'))]
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'filtered.maf'))]
            mtt.recordCommands([cmd], tmpDir, outPipes=outpipes)
            mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
            self.assertTrue(mafIsFiltered(os.path.join(tmpDir, 'filtered.maf'), expectedOutput))
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
            expectedOutput = []
            tmpDir = os.path.abspath(mtt.makeTempDir('memory1'))
            order = [1] * len(g_duplicateBlocks) + [0] * len(g_nonDuplicateBlocks)
            random.shuffle(order)
            random.shuffle(g_duplicateBlocks)
            random.shuffle(g_nonDuplicateBlocks)
            j, k = 0, 0
            for dupBlock in order:
                if dupBlock:
                    shuffledBlocks.append(g_duplicateBlocks[j][0])
                    expectedOutput.append(g_duplicateBlocks[j][1])
                    j += 1
                else:
                    shuffledBlocks.append(g_nonDuplicateBlocks[k])
                    expectedOutput.append(g_nonDuplicateBlocks[k])
                    k += 1
            testMaf = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test.maf')),
                                   ''.join(shuffledBlocks), g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = mtt.genericValgrind(tmpDir)
            cmd += [os.path.abspath(os.path.join(parent, 'test', 'mafDuplicateFilter')), 
                   '--maf', os.path.abspath(os.path.join(tmpDir, 'test.maf'))]
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'filtered.maf'))]
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
        for i in xrange(0, 10):
            tmpDir = os.path.abspath(mtt.makeTempDir('memory2'))
            random.shuffle(g_nonDuplicateBlocks)
            testMaf = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test.maf')),
                                   ''.join(g_nonDuplicateBlocks), g_headers)
            expectedOutput = g_nonDuplicateBlocks
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = mtt.genericValgrind(tmpDir)
            cmd += [os.path.abspath(os.path.join(parent, 'test', 'mafDuplicateFilter')), 
                   '--maf', os.path.abspath(os.path.join(tmpDir, 'test.maf'))]
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'filtered.maf'))]
            mtt.recordCommands([cmd], tmpDir, outPipes=outpipes)
            mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
            self.assertTrue(mtt.noMemoryErrors(os.path.join(tmpDir, 'valgrind.xml')))
            mtt.removeDir(tmpDir)

if __name__ == '__main__':
    unittest.main()
