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
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), '../../mafValidator/src/')))
import mafValidator as mafval

g_headers = ['''##maf version=1 scoring=tba.v8
# tba.v8 (((human chimp) baboon) (mouse rat))

''',
             '''##maf version=1 scoring=tba.v8
# tba.v8 (((human chimp) baboon) (mouse rat))
''']

g_blocks = [('''a score=0
#block 0, when target is in correct strand, do nothing
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
#block 0, when target is in correct strand, do nothing
s target.chr0        38 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 gcagctgaaaaca
s name.chr1           0 13 +       100 gcagctgaaaaca
s name2.chr1         50 13 +       100 gcagctgaaaaca
s name3.chr9         50 13 +       100 gcagctgaaaaca
s name4.chr&         50 13 +       100 gcagctgaaaaca
s name4.chrA         50 13 +       100 gcagctgaaaacT

''', '+'),
                     ('''a score=0
#block 1, when target is in reverse strand flip block
s name                0 13 +       100 gcagctgaaaaca
s name2.chr1         50 13 +       100 gcagctgaaaaca
s name3.chr9         50 13 +       100 gcagctgaaaaca
s name4.chr&         50 13 +       100 gcagctgaaaaca
s name5              50 13 +       100 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s target.chr0 158545457 13 - 158545518 gcagctgaaaacT
s target.chr1 158545457 13 + 158545518 gcagctgaaaaca

''',
                      '''a score=0
#block 1, when target is in reverse strand flip block
s name                87 13 -       100 tgttttcagctgc
s name2.chr1          37 13 -       100 tgttttcagctgc
s name3.chr9          37 13 -       100 tgttttcagctgc
s name4.chr&          37 13 -       100 tgttttcagctgc
s name5               37 13 -       100 tgttttcagctgc
s panTro1.chr6 132707175 13 - 161576975 tgttttcagctgc
s baboon         4373603 13 -   4622798 tgttttcagctgc
s target.chr0         48 13 + 158545518 Agttttcagctgc
s target.chr1         48 13 - 158545518 tgttttcagctgc

''', '+'),
                     ('''a score=0
#block 2, when the target is in correct strand, do nothing
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
#block 2, when the target is in correct strand, do nothing
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

''', '+'),
                     ('''a score=0
#block 3, when target is in reverse strand, flip block
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
#block 3, when target is in reverse strand, flip block
s name                77 13 -       100 tgttttcagctgc
s name.chr1           37 13 -       100 agttttcagctgc
s name.chr2           37 13 -       100 agttttcagctgc
s name3.chr9          37 13 -       100 tgttttcagctgc
s name4.chr&          37 13 -       100 tgttttcagctgc
s name5               37 13 -       100 tgttttcagctgc
s panTro1.chr6 132707175 13 - 161576975 agtattcagctgT
s target.chr0  158545447  9 - 158545518 tgt-ttc---tgc
s baboon         4373603 13 -   4622798 tgttttcagctgc
s hg16.chr7    130838284 13 - 158545518 tgttttcagctgc
s panTro1.chr7 132707175 13 - 161576975 Agtattcagctgc
s baboon         4373603 13 -   4622798 tgttttcagctgc
s mm4.chr6      97794610 13 - 151104725 TATTTTCAGCTGT

''', '-'),
                     ('''a score=0
#block 4, when the target is in correct strand, do nothing
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
#block 4, when the target is in correct strand, do nothing
s target.chr0        38 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 gcagctgannnnn
s name.chr1           0 13 +       100 gcagctgaaaacN
s name2.chr1         50 13 +       100 gcagctgaaaacN
s name3.chr9         50 13 +       100 gcagctgaaaaca
s name4.chr&         50 13 +       100 gcagctgaaaaca
s name4.chrA         50 13 +       100 gcagctgaaaacT

''', '+'),
            ('''a score=0
#block 5, when the target is present twice with both strands, do nothing
s name                0 13 +       100 gcagctgaaaaca
s name2.chr1         50 13 +       100 gcagctgaaaaca
s name3.chr9         50 13 +       100 gcagctgaaaaca
s name4.chr&         50 13 +       100 gcagctgaaaaca
s name5              50 13 +       100 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s target.chr0 158545457 13 - 158545518 gcagctgaaaacT
s target.chr0 158545407 13 + 158545518 gcagctgaaaaca

''',
             '''a score=0
#block 5, when the target is present twice with both strands, do nothing
s name                0 13 +       100 gcagctgaaaaca
s name2.chr1         50 13 +       100 gcagctgaaaaca
s name3.chr9         50 13 +       100 gcagctgaaaaca
s name4.chr&         50 13 +       100 gcagctgaaaaca
s name5              50 13 +       100 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s target.chr0 158545457 13 - 158545518 gcagctgaaaacT
s target.chr0 158545407 13 + 158545518 gcagctgaaaaca

''', '+'),
            ]
class GenericValidationOptions:
    def __init__(self):
        self.lookForDuplicateColumns = False
        self.testChromNames = False
        self.validateSequence = True
def hashify(s):
    return s.replace(' ', '').replace('\n', '')
def mafIsCoerced(maf, expected):
    f = open(maf)
    lastLine = mtt.processHeader(f)
    # walk through the maf, assessing the equivalence to the blockList items
    b = mtt.extractBlockStr(f, lastLine)
    lastLine = None
    if hashify(b) != hashify(expected):
        print 'dang'
        print 'observed:'
        print b
        print '!='
        print 'expected:'
        print expected
        return False
    return True
class CoercionTest(unittest.TestCase):
    def testCoercion(self):
        """ mafBlockStrandCoercer should coerce blocks into strandedness based on user input
        """
        mtt.makeTempDirParent()
        tmpDir = os.path.abspath(mtt.makeTempDir('coercion'))
        customOpts = GenericValidationOptions()
        for a, expected, strand in g_blocks:
            testMaf = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test.maf')),
                                   ''.join(a), g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = [os.path.abspath(os.path.join(parent, 'test', 'mafBlockStrandCoercer')), 
                   '--maf', os.path.abspath(os.path.join(tmpDir, 'test.maf')),
                   '--seq', 'target.chr0', '--strand', strand]
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'coerced.maf'))]
            mtt.recordCommands([cmd], tmpDir, outPipes=outpipes)
            mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
            self.assertTrue(mafIsCoerced(os.path.join(tmpDir, 'coerced.maf'), expected))
            self.assertTrue(mafval.validateMaf(os.path.join(tmpDir, 'coerced.maf'), customOpts))
        mtt.removeDir(tmpDir)
    def testMemory1(self):
        """ If valgrind is installed on the system, check for memory related errors (1).
        """
        mtt.makeTempDirParent()
        valgrind = mtt.which('valgrind')
        if valgrind is None:
            return
        tmpDir = os.path.abspath(mtt.makeTempDir('memeory1'))
        customOpts = GenericValidationOptions()
        for a, expected, strand in g_blocks:
            testMaf = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test.maf')),
                                   ''.join(a), g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = mtt.genericValgrind(tmpDir)
            cmd += [os.path.abspath(os.path.join(parent, 'test', 'mafBlockStrandCoercer')), 
                    '--maf', os.path.abspath(os.path.join(tmpDir, 'test.maf')),
                    '--seq', 'target', '--strand', strand]
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'coerced.maf'))]
            mtt.recordCommands([cmd], tmpDir, outPipes=outpipes)
            mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
            self.assertTrue(mtt.noMemoryErrors(os.path.join(tmpDir, 'valgrind.xml')))
            self.assertTrue(mafval.validateMaf(os.path.join(tmpDir, 'coerced.maf'), customOpts))
        mtt.removeDir(tmpDir)

if __name__ == '__main__':
    unittest.main()
