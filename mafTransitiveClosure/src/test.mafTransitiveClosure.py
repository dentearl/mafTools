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

g_headers = ['''##maf version=1 scoring=tba.v8
# tba.v8 (((human chimp) baboon) (mouse rat))

''',
             '''##maf version=1 scoring=tba.v8
# tba.v8 (((human chimp) baboon) (mouse rat))
''']

def mafIsClosed(maf, outList):
    f = open(maf)
    lastLine = mtt.processHeader(f)
    for i in xrange(0, len(outList)):
        block = mtt.extractBlockStr(f, lastLine)
        lastLine = None
        block = block.split('\n')
        block = block[1:] # throw away alignment line
        for line in block:
            line = line.strip()
            if line == '':
                continue
            if line not in outList:
                return False
    f.close()
    return True
    
class TransitiveClosureTest(unittest.TestCase):
    knownResults = [('''a score=23262.0     
s hg18.chr7      0 16 + 100 ATTGTCTCTTACGGTG
s panTro1.chr6   0 16 + 100 ATTGTCTCTTACGGTG

a score=23262.0     
s hg18.chr7      0 16 + 100 ATTGTCTCTTACGGTG
s baboon         0 16 + 100 GTTGTCTCTTATGGTG

a score=23262.0       
s hg18.chr7      0 16 + 100 ATTGTCTCTTACGGTG
s mm4.chr6       0 16 + 100 ATTGTCTCTCAGTGTG

a score=23262.0     
s hg18.chr7      0 16 + 100 ATTGTCTCTTACGGTG
s rn3.chr4       1 16 + 100 GTTGTCTCTCAATGTG
                   
''',
                     ['a degree=5',
                      's hg18.chr7      0 16 + 100 ATTGTCTCTTACGGTG',
                      's panTro1.chr6   0 16 + 100 ATTGTCTCTTACGGTG',
                      's baboon         0 16 + 100 GTTGTCTCTTATGGTG',
                      's mm4.chr6       0 16 + 100 ATTGTCTCTCAGTGTG',
                      's rn3.chr4       1 16 + 100 GTTGTCTCTCAATGTG']),
                    ('''a score=23262.0       
s hg18.chr7      0 10 - 10 AAAAAGGGGG
s mm4.chr6       1 10 - 11 GATTGTCTCC

a score=23262.0     
s hg18.chr7      4 6 + 10 CTTTTT
s rn3.chr4       4 6 - 12 GGGGGT

''',
                     ['a degree=2',
                      's hg18.chr7   0 4 + 10 CCCC',
                      's mm4.chr6    0 4 + 11 GGAG',
                      ''
                      'a degree=3',
                      's hg18.chr7   4 6 + 10 CTTTTT',
                      's mm4.chr6    4 6 + 11 ACAATC',
                      's rn3.chr4    4 6 - 12 GGGGGT',
                      ]),
                    ('''a score=23262.0     
s hg18.chr7      0 10 + 20 AATTGTCTCT
s panTro1.chr6   0 10 + 20 AATTGTCTCT

a score=23262.0     
s hg18.chr7      10 10 + 20 CCCCCTTTTT
s panTro1.chr6   10 10 + 20 GGCCTTAATT

a score=23262.0     
s hg18.chr7      5 10 - 20 GGGGGAGAGA
s baboon         5 10 - 20 GGTTGTCTCT

a score=23262.0     
s hg18.chr7      0 10 - 20 AAAAAGGGGG
s baboon         0 10 - 20 AATTAGGTTG

a score=23262.0     
s hg18.chr7      10 10 - 20 AGAGACAATT
s baboon         10 10 - 20 TCTCTCCGGG

a score=23262.0       
s hg18.chr7      0 20 - 20 AAAAAGGGGGAGAGACAATT
s mm4.chr6       0 20 - 20 GATTGTCTCTAAATTTAAAT

a score=23262.0     
s hg18.chr7      10 10 + 20 CCCCCTTTTT
s rn3.chr4       10 10 - 20 GGTTGTCTCT

a score=23262.0     
s hg18.chr7      0 10 + 20 AATTGTCTCT
s rn3.chr4       0 10 - 20 CCGGCGGTTG

                   
''',
                     ['a degree=5',
                      's hg18.chr7      0 20 + 20 AATTGTCTCTCCCCCTTTTT',
                      's panTro1.chr6   0 20 + 20 AATTGTCTCTGGCCTTAATT',
                      's baboon         0 20 + 20 CCCGGAGAGACAACCTAATT',
                      's mm4.chr6       0 20 + 20 ATTTAAATTTAGAGACAATC',
                      's rn3.chr4       0 20 - 20 CCGGCGGTTGGGTTGTCTCT',]),]
    def testKnownInOut(self):
        """ mafTransitiveClosure should compute the transitive closure of a maf built by pairwise alignment to a reference sequence.
        """
        mtt.makeTempDirParent()
        tmpDir = os.path.abspath(mtt.makeTempDir('knownInOut'))
        for inMaf, outList in self.knownResults:
            testMaf = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test.maf')), 
                                   inMaf, g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = [os.path.abspath(os.path.join(parent, 'test', 'mafTransitiveClosure')), 
                   '--maf', os.path.abspath(os.path.join(tmpDir, 'test.maf'))]
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'transitiveClosure.maf'))]
            mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
            passed = mafIsClosed(os.path.join(tmpDir, 'transitiveClosure.maf'), outList)
            self.assertTrue(passed)
        mtt.removeDir(tmpDir)
    def testMemory1(self):
        """ mafTransitiveClosure should be memory clean.
        """
        mtt.makeTempDirParent()
        valgrind = mtt.which('valgrind')
        if valgrind is None:
            return
        tmpDir = os.path.abspath(mtt.makeTempDir('memory1'))
        parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        cmd = mtt.genericValgrind(tmpDir)
        cmd += [os.path.abspath(os.path.join(parent, 'test', 'mafTransitiveClosure')), 
               '--maf', os.path.abspath(os.path.join(os.curdir, 'test.maf'))]
        outpipes = [os.path.abspath(os.path.join(tmpDir, 'transitiveClosure.maf'))]
        mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
        passed = mtt.noMemoryErrors(os.path.join(tmpDir, 'valgrind.xml'))
        self.assertTrue(passed)
        mtt.removeDir(tmpDir)
    def testMemory2(self):
        """ the CuTest tests should be memory clean.
        """
        valgrind = mtt.which('valgrind')
        if valgrind is None:
            return
        tmpDir = os.path.abspath(mtt.makeTempDir('memory2'))
        parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        cmd = mtt.genericValgrind(tmpDir)
        cmd += [os.path.abspath(os.path.join(parent, 'test', 'mafTransitiveClosure')), '--test']
        outpipes = [os.path.abspath(os.path.join(tmpDir, 'transitiveClosure.maf'))]
        mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
        passed = mtt.noMemoryErrors(os.path.join(tmpDir, 'valgrind.xml'))
        self.assertTrue(passed)
        mtt.removeDir(tmpDir)

if __name__ == '__main__':
    unittest.main()
