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

g_target = 'hg18.chr7'
g_headers = ['''##maf version=1 scoring=tba.v8
# tba.v8 (((human chimp) baboon) (mouse rat))

''',
           '''##maf version=1 scoring=tba.v8
# tba.v8 (((human chimp) baboon) (mouse rat))
''',]

# does not contain target hg18.chr7
g_nonTargetBlocks = ['''a score=23261.0
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon.chr0    249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA

''',
                     '''a score=23261.0
s banana.chr6   28869787 13 + 161576975 gcagctgaaaaca
s apple.chr0      249182 13 +   4622798 gcagctgaaaaca
s fork.chr6     53310102 13 + 151104725 ACAGCTGAAAATA

''',
                   '''a score=0
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
i panTro1.chr6 N 0 C 0
s baboon.chr0    249182 13 +   4622798 gcagctgaaaaca
i baboon.chr0   I 234 n 19

''',
                   '''a score=23262.0     
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon.chr0    116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG

''',
                   ]

# target is hg18.chr7, blocks MUST BE in correct order.
g_targetBlocks = ['''a score=23263.0
# sorted block 0
s hg18.chr7     7578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon.chr0    116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
                   
''',
                  '''a score=23263.0
# sorted block 1
s hg18.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon.chr0    116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
                   
''',
                '''a score=5062.0                    
# sorted block 2
s hg18.chr7    27699739 6 + 158545518 TAAAGA
s panTro1.chr6 28862317 6 + 161576975 TAAAGA
s baboon.chr0    241163 6 +   4622798 TAAAGA 
s mm4.chr6     53303881 6 + 151104725 TAAAGA
s rn3.chr4     81444246 6 + 187371129 taagga

''',
                '''a score=23264.0
# sorted block 3
s hg18.chr7    27707000 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon.chr0    249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA

''',
                '''a score=6636.0
# sorted block 4
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon.chr0    249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA
s hg18.chr7    27707221 13 + 158545518 gcagctgaaaaca

''',
                '''a score=0
# sorted block 5
s hg18.chr7              237249719 26 - 247249719 TTTTTGAAAAACAAACAACAAGTTGG
s panTro2.chrUn            9697231 26 +  58616431 TTTTTGAAAAACAAACAACAAGTTGG
q panTro2.chrUn                                   99999999999999999999999999
s dasNov1.scaffold_179265     1474  7 +      4584 TT----------AAGCA---------
q dasNov1.scaffold_179265                         99----------32239--------- 

''',
                ]

def mafIsSorted(maf):
    f = open(maf)
    lastLine = mtt.processHeader(f)
    observedNonTargetBlocks = []
    expectedNonTargetBlocks = list(g_nonTargetBlocks)
    for i in xrange(0, len(expectedNonTargetBlocks)):
        observedNonTargetBlocks.append(mtt.extractBlockStr(f, lastLine).replace(' ', ''))
        lastLine = None
        expectedNonTargetBlocks[i] = expectedNonTargetBlocks[i].replace(' ', '')
    if observedNonTargetBlocks != expectedNonTargetBlocks:
        f.close()
        print '\nnon-target block failure observed:'
        print ''.join(observedNonTargetBlocks)
        print '!= expected:'
        print ''.join(expectedNonTargetBlocks)
        return False
    sortedBlocks = []
    expectedTargetBlocks = list(g_targetBlocks)
    for i in xrange(0, len(expectedTargetBlocks)):
        sortedBlocks.append(mtt.extractBlockStr(f).replace(' ', ''))
        expectedTargetBlocks[i] = expectedTargetBlocks[i].replace(' ', '')
    for i in xrange(0, len(sortedBlocks)):
        if sortedBlocks[i] != expectedTargetBlocks[i]:
            f.close()
            print '\nsorted block failure'
            print sortedBlocks[i]
            print '!='
            print expectedTargetBlocks[i]
            print 'observed '
            print sortedBlocks
            print '!= expected '
            print expectedTargetBlocks
            return False
    f.close()
    return True

class SortTest(unittest.TestCase):
    def testSorting(self):
        """ Blocks should be sorted by the start field of the target sequence, blocks that do not contain the target sequence should appear in the output at the start of the file, in the same order they appear in the input.
        """
        shuffledTargets = list(g_targetBlocks)
        for i in xrange(0, 200):
            tmpDir = os.path.abspath(mtt.makeTempDir('sorting'))
            random.shuffle(g_nonTargetBlocks)
            random.shuffle(shuffledTargets)
            shuffledBlocks = list(shuffledTargets)
            lower = 0
            for j in xrange(0, len(g_nonTargetBlocks)):
                # randomly insert the non target blocks, but keep a record
                # of their relative order.
                index = random.randint(lower, len(shuffledBlocks))
                shuffledBlocks.insert(index, g_nonTargetBlocks[j])
                lower = index + 1
            testMaf = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test.maf')), 
                                   ''.join(shuffledBlocks), g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = [os.path.abspath(os.path.join(parent, 'test', 'mafBlockSorter'))]
            cmd += ['--maf', os.path.abspath(os.path.join(tmpDir, 'test.maf')), 
                    '--seq', 'hg18.chr7']
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'sorted.maf'))]
            mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
            self.assertTrue(mafIsSorted(os.path.join(tmpDir, 'sorted.maf')))
            mtt.removeDir(tmpDir)
    def testMemory1(self):
        """ If valgrind is installed on the system, check for memory related errors (1).
        """
        valgrind = mtt.which('valgrind')
        if valgrind is None:
            return
        shuffledTargets = list(g_targetBlocks)
        for i in xrange(0, 20):
            tmpDir = os.path.abspath(mtt.makeTempDir('memory1'))
            random.shuffle(g_nonTargetBlocks)
            random.shuffle(shuffledTargets)
            shuffledBlocks = list(shuffledTargets)
            lower = 0
            for j in xrange(0, len(g_nonTargetBlocks)):
                # randomly insert the non target blocks, but keep a record
                # of their relative order.
                index = random.randint(lower, len(shuffledBlocks))
                shuffledBlocks.insert(index, g_nonTargetBlocks[j])
                lower = index + 1
            testMaf = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test.maf')),
                                   ''.join(shuffledBlocks), g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = mtt.genericValgrind(tmpDir)
            cmd.append(os.path.abspath(os.path.join(parent, 'test', 'mafBlockSorter')))
            cmd += ['--maf', os.path.abspath(os.path.join(tmpDir, 'test.maf')), 
                    '--seq', 'hg18.chr7']
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'sorted.maf'))]
            mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
            self.assertTrue(mtt.noMemoryErrors(os.path.join(tmpDir, 'valgrind.xml')))
            mtt.removeDir(tmpDir)

if __name__ == '__main__':
    unittest.main()
