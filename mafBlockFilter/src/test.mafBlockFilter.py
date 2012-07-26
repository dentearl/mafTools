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

g_sequenceList = 'target0,target1,target2'
g_header = None
g_headers = ['''##maf version=1 scoring=tba.v8
# tba.v8 (((human chimp) baboon) (mouse rat))

''', 
             '''##maf version=1 scoring=tba.v8
# tba.v8 (((human chimp) baboon) (mouse rat))
''',]

# these are triples of blocks, target positions for True, 
# line numbers where the target occurs.
g_knownIncludes = [('''a score=0
# test1
s target0.chr0         0 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA
s name.chr1           0 10 +       100 ATGT---ATGCCG
s name2.chr1         50 10 +       100 ATGT---ATGCCG
s name3.chr9         50 10 +       100 ATGTA---TGCCG
s name4.chr&         50 10 +       100 ATG---TATGCCG
s name5 50 10 + 100 ATGTATGCCG

''',
                    '''a score=0
# test1
s target0.chr0         0 13 + 158545518 gcagctgaaaaca

'''),
                   ('''a score=0
# test2
s target1.chr0        38 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA
s name.chr1           0 10 +       100 ATGT---ATGCCG
s name2.chr1         50 10 +       100 ATGT---ATGCCG
s name3.chr9         50 10 +       100 ATGTA---TGCCG
s name4.chr&         50 10 +       100 ATG---TATGCCG
s name5 50 10 + 100 ATGTATGCCG

''', 
                    '''a score=0
# test2
s target1.chr0        38 13 + 158545518 gcagctgaaaaca

''' ),
                   ('''a score=0
# test3
s name                0 10 +       100 ATGTATGC---CG
s name2.chr1         50 10 +       100 ATGTATG---CCG
s name3.chr9         50 10 +       100 ATGTATG---CCG
s name4.chr&         50 10 +       100 ATGTATG---CCG
s name5              50 10 +       100 ATGTATG---CCG
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s target2.chr0 158545457 10 - 158545518 ATGTATG---CCG

a score=0
s name                0 10 +       100 ATGTATGC---CG
s name2.chr1         50 10 +       100 ATGTATG---CCG
s name3.chr9         50 10 +       100 ATGTATG---CCG
s name4.chr&         50 10 +       100 ATGTATG---CCG

''', 
                    '''a score=0
# test3
s target2.chr0 158545457 10 - 158545518 ATGTATG---CCG

'''),

                   ('''a score=0
# test4
s name               10 10 +       100 ATGTAT---GCCG
s name2.chr1         50 10 +       100 ATGTAT---GCCG
s name3.chr9         50 10 +       100 ATGTAT---GCCG
s name4.chr&         50 10 +       100 ATGTAT---GCCG
s name5              50 10 +       100 ATGTAT---GCCG
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s target0.chr0        62  9 + 158545518 gca---gaa-aca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s hg16.chr7    27707221 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA

a score=0
s name                0 10 +       100 ATGTATGC---CG
s name2.chr1         50 10 +       100 ATGTATG---CCG
s name3.chr9         50 10 +       100 ATGTATG---CCG
s name4.chr&         50 10 +       100 ATGTATG---CCG

a score=0
s name                0 10 +       100 ATGTATGC---CG
s name2.chr1         50 10 +       100 ATGTATG---CCG
s name3.chr9         50 10 +       100 ATGTATG---CCG
s name4.chr&         50 10 +       100 ATGTATG---CCG

a score=0
s name                0 10 +       100 ATGTATGC---CG
s target1.chr0       50 10 +       100 ATGTATG---CCG
s name2.chr1         50 10 +       100 ATGTATG---CCG
s name3.chr9         50 10 +       100 ATGTATG---CCG
s name4.chr&         50 10 +       100 ATGTATG---CCG
''', 
                    '''a score=0
# test4
s target0.chr0        62  9 + 158545518 gca---gaa-aca

a score=0
s target1.chr0       50 10 +       100 ATGTATG---CCG

'''),
            # the overlap is 9 bases long, 62-70
            ]
g_knownExcludes = [('''a score=0
# test5
s target.chr0         0 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA
s name.chr1           0 10 +       100 ATGT---ATGCCG
s name2.chr1         50 10 +       100 ATGT---ATGCCG
s name3.chr9         50 10 +       100 ATGTA---TGCCG
s name4.chr&         50 10 +       100 ATG---TATGCCG
s name5 50 10 + 100 ATGTATGCCG

''',
                    '''a score=0
# test5
s target.chr0         0 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA
s name.chr1           0 10 +       100 ATGT---ATGCCG
s name2.chr1         50 10 +       100 ATGT---ATGCCG
s name3.chr9         50 10 +       100 ATGTA---TGCCG
s name4.chr&         50 10 +       100 ATG---TATGCCG
s name5 50 10 + 100 ATGTATGCCG

'''),
                   ('''a score=0
# test6
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA
s name.chr1           0 10 +       100 ATGT---ATGCCG
s name2.chr1         50 10 +       100 ATGT---ATGCCG
s name3.chr9         50 10 +       100 ATGTA---TGCCG
s name4.chr&         50 10 +       100 ATG---TATGCCG
s name5 50 10 + 100 ATGTATGCCG

''', 
                    '''a score=0
# test6
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA
s name.chr1           0 10 +       100 ATGT---ATGCCG
s name2.chr1         50 10 +       100 ATGT---ATGCCG
s name3.chr9         50 10 +       100 ATGTA---TGCCG
s name4.chr&         50 10 +       100 ATG---TATGCCG
s name5 50 10 + 100 ATGTATGCCG

''' ),
                   ('''a score=0
# test7
s name                0 10 +       100 ATGTATGC---CG
s name2.chr1         50 10 +       100 ATGTATG---CCG
s name3.chr9         50 10 +       100 ATGTATG---CCG
s name4.chr&         50 10 +       100 ATGTATG---CCG
s name5              50 10 +       100 ATGTATG---CCG
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s target2.chr0 158545457 10 - 158545518 ATGTATG---CCG

a score=0
s name                0 10 +       100 ATGTATGC---CG
s name2.chr1         50 10 +       100 ATGTATG---CCG
s name3.chr9         50 10 +       100 ATGTATG---CCG
s name4.chr&         50 10 +       100 ATGTATG---CCG

''', 
                    '''a score=0
# test7
s name                0 10 +       100 ATGTATGC---CG
s name2.chr1         50 10 +       100 ATGTATG---CCG
s name3.chr9         50 10 +       100 ATGTATG---CCG
s name4.chr&         50 10 +       100 ATGTATG---CCG
s name5              50 10 +       100 ATGTATG---CCG
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca

a score=0
s name                0 10 +       100 ATGTATGC---CG
s name2.chr1         50 10 +       100 ATGTATG---CCG
s name3.chr9         50 10 +       100 ATGTATG---CCG
s name4.chr&         50 10 +       100 ATGTATG---CCG

'''),

                   ('''a score=0
# test8
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

a score=0
s name                0 10 +       100 ATGTATGC---CG
s name2.chr1         50 10 +       100 ATGTATG---CCG
s name3.chr9         50 10 +       100 ATGTATG---CCG
s name4.chr&         50 10 +       100 ATGTATG---CCG

a score=0
s name                0 10 +       100 ATGTATGC---CG
s name2.chr1         50 10 +       100 ATGTATG---CCG
s name3.chr9         50 10 +       100 ATGTATG---CCG
s name4.chr&         50 10 +       100 ATGTATG---CCG

a score=0
s name                0 10 +       100 ATGTATGC---CG
s target2.chr0       50 10 +       100 ATGTATG---CCG
s name2.chr1         50 10 +       100 ATGTATG---CCG
s name3.chr9         50 10 +       100 ATGTATG---CCG
s name4.chr&         50 10 +       100 ATGTATG---CCG

''', 
                    '''a score=0
# test8
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

a score=0
s name                0 10 +       100 ATGTATGC---CG
s name2.chr1         50 10 +       100 ATGTATG---CCG
s name3.chr9         50 10 +       100 ATGTATG---CCG
s name4.chr&         50 10 +       100 ATGTATG---CCG

a score=0
s name                0 10 +       100 ATGTATGC---CG
s name2.chr1         50 10 +       100 ATGTATG---CCG
s name3.chr9         50 10 +       100 ATGTATG---CCG
s name4.chr&         50 10 +       100 ATGTATG---CCG

a score=0
s name                0 10 +       100 ATGTATGC---CG
s name2.chr1         50 10 +       100 ATGTATG---CCG
s name3.chr9         50 10 +       100 ATGTATG---CCG
s name4.chr&         50 10 +       100 ATGTATG---CCG

'''),
            # the overlap is 9 bases long, 62-70
            ]
g_knownDegreeLT = [('''a score=0
# test0
s target.chr0         0 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA
s name.chr1           0 10 +       100 ATGT---ATGCCG
s name2.chr1         50 10 +       100 ATGT---ATGCCG
s name3.chr9         50 10 +       100 ATGTA---TGCCG
s name4.chr&         50 10 +       100 ATG---TATGCCG
s name5 50 10 + 100 ATGTATGCCG

a score=0
# test0
s target.chr0         0 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s name5 50 10 + 100 ATGTATGCCG

a score=0
# test0
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA
s name5 50 10 + 100 ATGTATGCCG

''',
                    '''a score=0
# test0
s target.chr0         0 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA
s name.chr1           0 10 +       100 ATGT---ATGCCG
s name2.chr1         50 10 +       100 ATGT---ATGCCG
s name3.chr9         50 10 +       100 ATGTA---TGCCG
s name4.chr&         50 10 +       100 ATG---TATGCCG
s name5 50 10 + 100 ATGTATGCCG

a score=0
# test0
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA
s name5 50 10 + 100 ATGTATGCCG

'''),]

g_knownDegreeGT = [('''a score=0
# test0
s target.chr0         0 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA
s name.chr1           0 10 +       100 ATGT---ATGCCG
s name2.chr1         50 10 +       100 ATGT---ATGCCG
s name3.chr9         50 10 +       100 ATGTA---TGCCG
s name4.chr&         50 10 +       100 ATG---TATGCCG
s name5 50 10 + 100 ATGTATGCCG

a score=0
# test0
s target.chr0         0 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s name5 50 10 + 100 ATGTATGCCG

a score=0
# test0
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA
s name5 50 10 + 100 ATGTATGCCG

''',
                    '''a score=0
# test0
s target.chr0         0 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s name5 50 10 + 100 ATGTATGCCG

a score=0
# test0
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA
s name5 50 10 + 100 ATGTATGCCG

'''),]

def mafIsFiltered(filename, expected, header):
    f = open(filename)
    lastLine = mtt.processHeader(f)
    maf = ''
    b = mtt.extractBlockStr(f, lastLine)
    while (b is not None):
        lastLine = None
        maf += b
        b = mtt.extractBlockStr(f, lastLine)
    if maf != expected:
        print 'dang'
        print 'observed:'
        print maf
        print '!='
        print 'expected:'
        print expected
        f.close()
        return False
    f.close()
    return True

class FilterTest(unittest.TestCase):
    def testFilterIncludes(self):
        """ mafBlockFilter should report blocks that match the filter settings for --includeSeq.
        """
        global g_header
        mtt.makeTempDirParent()
        for i in xrange(0, len(g_knownIncludes)):
            tmpDir = os.path.abspath(mtt.makeTempDir('filterIncludes'))
            testMafPath, g_header = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test.maf')),
                                                 g_knownIncludes[i][0], g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = []
            cmd.append(os.path.abspath(os.path.join(parent, 'test', 'mafBlockFilter')))
            cmd += ['--maf', testMafPath, '--includeSeq', '%s' % g_sequenceList]
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'filtered.maf'))]
            mtt.recordCommands([cmd], tmpDir, outPipes=outpipes)
            mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
            filtered = mafIsFiltered(os.path.join(tmpDir, 'filtered.maf'), g_knownIncludes[i][1], g_header)
            self.assertTrue(filtered)
            if filtered:
                mtt.removeDir(tmpDir)
    def testFilterExcludes(self):
        """ mafBlockFilter should report blocks that match the filter settings for --excludeSeq.
        """
        global g_header
        mtt.makeTempDirParent()
        for i in xrange(0, len(g_knownExcludes)):
            tmpDir = os.path.abspath(mtt.makeTempDir('filterExcludes'))
            testMafPath, g_header = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test.maf')),
                                                 g_knownExcludes[i][0], g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = []
            cmd.append(os.path.abspath(os.path.join(parent, 'test', 'mafBlockFilter')))
            cmd += ['--maf', testMafPath, '--excludeSeq', '%s' % g_sequenceList]
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'filtered.maf'))]
            mtt.recordCommands([cmd], tmpDir, outPipes=outpipes)
            mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
            filtered = mafIsFiltered(os.path.join(tmpDir, 'filtered.maf'), g_knownExcludes[i][1], g_header)
            self.assertTrue(filtered)
            if filtered:
                mtt.removeDir(tmpDir)
    def testFilterDegreeLT(self):
        """ mafBlockFilter should report blocks that match the filter settings for --noDegreeLT.
        """
        global g_header
        mtt.makeTempDirParent()
        for i in xrange(0, len(g_knownDegreeLT)):
            tmpDir = os.path.abspath(mtt.makeTempDir('filterDegreeLT'))
            testMafPath, g_header = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test.maf')),
                                                 g_knownDegreeLT[i][0], g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = []
            cmd.append(os.path.abspath(os.path.join(parent, 'test', 'mafBlockFilter')))
            cmd += ['--maf', testMafPath, '--noDegreeLT=4']
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'filtered.maf'))]
            mtt.recordCommands([cmd], tmpDir, outPipes=outpipes)
            mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
            filtered = mafIsFiltered(os.path.join(tmpDir, 'filtered.maf'), g_knownDegreeLT[i][1], g_header)
            self.assertTrue(filtered)
            if filtered:
                mtt.removeDir(tmpDir)
    def testFilterDegreeGT(self):
        """ mafBlockFilter should report blocks that match the filter settings for --noDegreeGT.
        """
        global g_header
        mtt.makeTempDirParent()
        for i in xrange(0, len(g_knownDegreeGT)):
            tmpDir = os.path.abspath(mtt.makeTempDir('filterDegreeGT'))
            testMafPath, g_header = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test.maf')),
                                                 g_knownDegreeGT[i][0], g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = []
            cmd.append(os.path.abspath(os.path.join(parent, 'test', 'mafBlockFilter')))
            cmd += ['--maf', testMafPath, '--noDegreeGT=4']
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'filtered.maf'))]
            mtt.recordCommands([cmd], tmpDir, outPipes=outpipes)
            mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
            filtered = mafIsFiltered(os.path.join(tmpDir, 'filtered.maf'), g_knownDegreeGT[i][1], g_header)
            self.assertTrue(filtered)
            if filtered:
                mtt.removeDir(tmpDir)
    def testMemory1(self):
        """ If valgrind is installed on the system, check for memory related errors (1).
        """
        mtt.makeTempDirParent()
        valgrind = mtt.which('valgrind')
        if valgrind is None:
            return
        global g_header
        for i in xrange(0, len(g_knownIncludes)):
            tmpDir = os.path.abspath(mtt.makeTempDir('memory1'))
            testMafPath, g_header = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test.maf')),
                                                 g_knownIncludes[i][0], g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = mtt.genericValgrind(tmpDir)
            cmd.append(os.path.abspath(os.path.join(parent, 'test', 'mafBlockFilter')))
            cmd += ['--maf', testMafPath, '--includeSeq', '%s' % g_sequenceList]
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
        global g_header
        for i in xrange(0, len(g_knownExcludes)):
            tmpDir = os.path.abspath(mtt.makeTempDir('memory2'))
            testMafPath, g_header = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test.maf')),
                                                 g_knownExcludes[i][0], g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = mtt.genericValgrind(tmpDir)
            cmd.append(os.path.abspath(os.path.join(parent, 'test', 'mafBlockFilter')))
            cmd += ['--maf', testMafPath, '--excludeSeq', '%s' % g_sequenceList]
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'filtered.maf'))]
            mtt.recordCommands([cmd], tmpDir, outPipes=outpipes)
            mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
            self.assertTrue(mtt.noMemoryErrors(os.path.join(tmpDir, 'valgrind.xml')))
            mtt.removeDir(tmpDir)
    def testMemory3(self):
        """ If valgrind is installed on the system, check for memory related errors (2).
        """
        mtt.makeTempDirParent()
        valgrind = mtt.which('valgrind')
        if valgrind is None:
            return
        global g_header
        mtt.makeTempDirParent()
        for i in xrange(0, len(g_knownDegreeLT)):
            tmpDir = os.path.abspath(mtt.makeTempDir('memory3'))
            testMafPath, g_header = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test.maf')),
                                                 g_knownDegreeLT[i][0], g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = mtt.genericValgrind(tmpDir)
            cmd.append(os.path.abspath(os.path.join(parent, 'test', 'mafBlockFilter')))
            cmd += ['--maf', testMafPath, '--noDegreeLT=4']
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'filtered.maf'))]
            mtt.recordCommands([cmd], tmpDir, outPipes=outpipes)
            mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
            self.assertTrue(mtt.noMemoryErrors(os.path.join(tmpDir, 'valgrind.xml')))
            mtt.removeDir(tmpDir)
    def testMemory4(self):
        """ If valgrind is installed on the system, check for memory related errors (2).
        """
        mtt.makeTempDirParent()
        valgrind = mtt.which('valgrind')
        if valgrind is None:
            return
        global g_header
        mtt.makeTempDirParent()
        for i in xrange(0, len(g_knownDegreeGT)):
            tmpDir = os.path.abspath(mtt.makeTempDir('memory4'))
            testMafPath, g_header = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test.maf')),
                                                 g_knownDegreeGT[i][0], g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = mtt.genericValgrind(tmpDir)
            cmd.append(os.path.abspath(os.path.join(parent, 'test', 'mafBlockFilter')))
            cmd += ['--maf', testMafPath, '--noDegreeGT=4']
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'filtered.maf'))]
            mtt.recordCommands([cmd], tmpDir, outPipes=outpipes)
            mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
            self.assertTrue(mtt.noMemoryErrors(os.path.join(tmpDir, 'valgrind.xml')))
            mtt.removeDir(tmpDir)
if __name__ == '__main__':
    unittest.main()
