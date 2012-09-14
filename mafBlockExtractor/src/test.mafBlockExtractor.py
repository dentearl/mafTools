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
s name5              50 10 +       100 ATGT---ATGCCG

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
                          '''a degree=4
s simChimp.chrA     15393341 348 + 53121445 ATATTGAGGAGCAGGATGGGTATAGAGGCCCTGACCCATTAATGTGTAAGCACTAGGCAGCTGGGAGATACCCCAGAGGGCGGGGTCACTGAATTCACTGGCCCACCACTGTAAATACATTCTAACCAGTGGGTTTAGGGCTCTGTGCATTAGAACCACTCTGAAGAAGTGTAACACACCACCTAGTGAGCTGCCGGGCCGCCAGCAACTTCTTTTTCCCACATGACCCATGCAAGCCCGTGATTTCTCCCTGGTACATGATATTTGGGATTCCAGGGACCTAATGGAGCATGCTATTCCTGTGTTAGTTATCACTTCGAAGGGGGTGCAAGAGTGTAAGTAATGGGT
s simGorilla.chrA   15595743 348 + 53120926 ATATTGAGGAGCAGGATGGGTATAGAAGCCCTGACCCATTAATGTGTAAGCACTAGGCAGCTGGGAGATACTCCAGAGGGAGGGGTCACTGAATTCACTGGCCCACCACTGTAAATACATTCTAACCAGTGGGTTTAGGGCTCGGTGCATTAGAACCACCCTGAAGAAGTGTAACGCACCACCTAGTGAGCTGCCGGGCCGCCAGCAACTTCTTTTTCCCACATGACCCATGCATGCCCGTGATTTCTCCCTGGTACATGGTTTTTGGGATTCCAGGGACCTAATGGAGCATACTATTCCTGTGTTAGTTATCACTTCGAAGGGGGTGCGAGAGTGTAAGTAATGGGT
s simHuman.chrA     36713600 348 - 53106993 ATATTGAGGAGCAGGATGGGTATAGAAGCCCTGACCTATTAATGTGTAAGCACTAGGCAGCTGGGCGATACCCCAGAGGGAGGGGTCACTGAATTCACTGGCCCACCACTGTAAATACATTCTAACCAGTGGGTTTAGGGCTCTGTGCATTAGAACCACCCTGAAGAAGAGTAACGCACCACCTAGTGAGCTGCCGGGCCGCCAGCAAGTTCTTTTTCCCACATGACCCATGCAAGCCCGTGATTTCTCCCTGGTACATGATATTTGGGATTCCAGGGACCTAATGGAGCATGCTATTCCTGTGTTAGTTATCACTTCGAAGGGGGTGCAAGAGTGTAAGTAATGGGT
s simOrang.chrE       126390 348 + 37692687 ATTTTGAGGAGCAGGATGGGTATAGAAGCCCTGACCCATTAATGTGTGAGCTCTAGGCAGCTTGGAGATACTGCAGAGGGAGGGGTCACTGAATTCACTGGCCCACCACTGTAAATACATACTAACCGGTGGGTTTAGGGCTCTGTGCATTAGAACCACCCTGAGGAAGTGTAACGCACCACCTAGTGAGCTGCCGGGCCACCAGCAACTTCTTTTTCCCACATGACCCATGCAAGCCCGTGATTTCTCCCTGGTACATGATCTTTGGGATTCCAGGGACCTAATGGCGGATGCTATTCCTGTGTTAGTTATCACTTCGAAGGGGGCGCAAGAGTGTAAGTAATGGGT
s target.chr0     36713600 348 - 53106993 ATATTGAGGAGCAGGATGGGTATAGAAGCCCTGACCTATTAATGTGTAAGCACTAGGCAGCTGGGCGATACCCCAGAGGGAGGGGTCACTGAATTCACTGGCCCACCACTGTAAATACATTCTAACCAGTGGGTTTAGGGCTCTGTGCATTAGAACCACCCTGAAGAAGAGTAACGCACCACCTAGTGAGCTGCCGGGCCGCCAGCAAGTTCTTTTTCCCACATGACCCATGCAAGCCCGTGATTTCTCCCTGGTACATGATATTTGGGATTCCAGGGACCTAATGGAGCATGCTATTCCTGTGTTAGTTATCACTTCGAAGGGGGTGCAAGAGTGTAAGTAATGGGT

'''
    ]
class GenericValidationOptions:
    def __init__(self):
        self.lookForDuplicateColumns = False
        self.testChromNames = False
        self.validateSequence = True
def hashStr(s):
    return s.replace(' ', '').replace('\n', '')
g_overlappingBlocksHashed = set()
for b in g_overlappingBlocks:
    g_overlappingBlocksHashed.add(hashStr(b))
def mafIsExtracted(maf):
    f = open(maf)
    lastLine = mtt.processHeader(f)
    for i in xrange(0, len(g_overlappingBlocks)):
        b = mtt.extractBlockStr(f, lastLine)
        lastLine = None
        if hashStr(b) not in g_overlappingBlocksHashed:
            print 'dang, block'
            print b
            print 'is not in set of expected output'
            print g_overlappingBlocks
            f.close()
            return False
    f.close()
    return True
def mafIsHardExtracted(name, outblocks, maf):
    blockSet = set()
    for b in outblocks:
        blockSet.add(b.hashify())
    f = open(maf)
    lastLine = mtt.processHeader(f)
    r = mtt.extractBlockStr(f, lastLine)
    observedBlocks = set()
    while r is not None:
        b = rawBlockToObj(r)
        if b is not None:
            if b.hashify() not in blockSet:
                print '\n[%s]' % name
                print 'dang, hashed output contains'
                print b.hashify()
                print 'not in hashed expected output'
                print blockSet
                f.close()
                return False
            else:
                observedBlocks.add(b.hashify())
        r = mtt.extractBlockStr(f)
    for b in blockSet:
        if b not in observedBlocks:
            print '\n[%s]' % name
            print 'dang, expected hashed output block'
            print b
            print 'not observed in hashed output blocks'
            for block in observedBlocks:
                print block
            f.close()
            return False
    f.close()        
    return True

class ExtractionTest(unittest.TestCase):
    def testExtraction(self):
        """ mafBlockExtractor should output blocks that meet the criteria for extraction. That is they contain the target sequence and have at least one base in the target range.
        """
        mtt.makeTempDirParent()
        customOpts = GenericValidationOptions()
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
            testMaf = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test_%d.maf' % i)),
                                   ''.join(shuffledBlocks), g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = [os.path.abspath(os.path.join(parent, 'test', 'mafBlockExtractor'))]
            cmd += ['--maf', os.path.abspath(os.path.join(tmpDir, 'test_%d.maf' % i)),
                    '--seq', g_targetSeq, '--start', '%d' % g_targetRange[0], 
                    '--stop', '%d' % g_targetRange[1],
                    '--soft'
                    ]
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'extracted.maf'))]
            mtt.recordCommands([cmd], tmpDir, outPipes=outpipes)
            mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
            self.assertTrue(mafIsExtracted(os.path.join(tmpDir, 'extracted.maf')))
            self.assertTrue(mafval.validateMaf(os.path.join(tmpDir, 'extracted.maf'), customOpts))
            mtt.removeDir(tmpDir)
    def testNonExtraction0(self):
        """ mafBlockExtractor should not extract blocks when they do not match.
        """
        customOpts = GenericValidationOptions()
        mtt.makeTempDirParent()
        for i in xrange(0, 10):
            tmpDir = os.path.abspath(mtt.makeTempDir('nonExtraction_0'))
            random.shuffle(g_nonOverlappingBlocks)
            testMaf = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test_%d.maf' % i)),
                                   ''.join(g_nonOverlappingBlocks), g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = [os.path.abspath(os.path.join(parent, 'test', 'mafBlockExtractor'))]
            cmd += ['--maf', os.path.abspath(os.path.join(tmpDir, 'test_%d.maf' % i)),
                    '--seq', g_targetSeq, '--start', '%d' % g_targetRange[0], 
                    '--stop', '%d' % g_targetRange[1],
                    '--soft'
                    ]
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'extracted.maf'))]
            mtt.recordCommands([cmd], tmpDir, outPipes=outpipes)
            mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
            self.assertTrue(mtt.mafIsEmpty(os.path.join(tmpDir, 'extracted.maf'), g_headers))
            self.assertTrue(mafval.validateMaf(os.path.join(tmpDir, 'extracted.maf'), customOpts))
            mtt.removeDir(tmpDir)
    def testNonExtraction1(self):
        """ mafBlockExtractor should not extract blocks when they do not match.
        """
        customOpts = GenericValidationOptions()
        mtt.makeTempDirParent()
        tmpDir = os.path.abspath(mtt.makeTempDir('nonExtraction_1'))
        block = '''a degree=4
s simChimp.chrA     15393341 348 + 53121445 ATATTGAGGAGCAGGATGGGTATAGAGGCCCTGACCCATTAATGTGTAAGCACTAGGCAGCTGGGAGATACCCCAGAGGGCGGGGTCACTGAATTCACTGGCCCACCACTGTAAATACATTCTAACCAGTGGGTTTAGGGCTCTGTGCATTAGAACCACTCTGAAGAAGTGTAACACACCACCTAGTGAGCTGCCGGGCCGCCAGCAACTTCTTTTTCCCACATGACCCATGCAAGCCCGTGATTTCTCCCTGGTACATGATATTTGGGATTCCAGGGACCTAATGGAGCATGCTATTCCTGTGTTAGTTATCACTTCGAAGGGGGTGCAAGAGTGTAAGTAATGGGT
s simGorilla.chrA   15595743 348 + 53120926 ATATTGAGGAGCAGGATGGGTATAGAAGCCCTGACCCATTAATGTGTAAGCACTAGGCAGCTGGGAGATACTCCAGAGGGAGGGGTCACTGAATTCACTGGCCCACCACTGTAAATACATTCTAACCAGTGGGTTTAGGGCTCGGTGCATTAGAACCACCCTGAAGAAGTGTAACGCACCACCTAGTGAGCTGCCGGGCCGCCAGCAACTTCTTTTTCCCACATGACCCATGCATGCCCGTGATTTCTCCCTGGTACATGGTTTTTGGGATTCCAGGGACCTAATGGAGCATACTATTCCTGTGTTAGTTATCACTTCGAAGGGGGTGCGAGAGTGTAAGTAATGGGT
s simHuman.chrA     36713600 348 - 53106993 ATATTGAGGAGCAGGATGGGTATAGAAGCCCTGACCTATTAATGTGTAAGCACTAGGCAGCTGGGCGATACCCCAGAGGGAGGGGTCACTGAATTCACTGGCCCACCACTGTAAATACATTCTAACCAGTGGGTTTAGGGCTCTGTGCATTAGAACCACCCTGAAGAAGAGTAACGCACCACCTAGTGAGCTGCCGGGCCGCCAGCAAGTTCTTTTTCCCACATGACCCATGCAAGCCCGTGATTTCTCCCTGGTACATGATATTTGGGATTCCAGGGACCTAATGGAGCATGCTATTCCTGTGTTAGTTATCACTTCGAAGGGGGTGCAAGAGTGTAAGTAATGGGT
s simOrang.chrE       126390 348 + 37692687 ATTTTGAGGAGCAGGATGGGTATAGAAGCCCTGACCCATTAATGTGTGAGCTCTAGGCAGCTTGGAGATACTGCAGAGGGAGGGGTCACTGAATTCACTGGCCCACCACTGTAAATACATACTAACCGGTGGGTTTAGGGCTCTGTGCATTAGAACCACCCTGAGGAAGTGTAACGCACCACCTAGTGAGCTGCCGGGCCACCAGCAACTTCTTTTTCCCACATGACCCATGCAAGCCCGTGATTTCTCCCTGGTACATGATCTTTGGGATTCCAGGGACCTAATGGCGGATGCTATTCCTGTGTTAGTTATCACTTCGAAGGGGGCGCAAGAGTGTAAGTAATGGGT
s target.chr0       36713600 348 - 53106993 ATATTGAGGAGCAGGATGGGTATAGAAGCCCTGACCTATTAATGTGTAAGCACTAGGCAGCTGGGCGATACCCCAGAGGGAGGGGTCACTGAATTCACTGGCCCACCACTGTAAATACATTCTAACCAGTGGGTTTAGGGCTCTGTGCATTAGAACCACCCTGAAGAAGAGTAACGCACCACCTAGTGAGCTGCCGGGCCGCCAGCAAGTTCTTTTTCCCACATGACCCATGCAAGCCCGTGATTTCTCCCTGGTACATGATATTTGGGATTCCAGGGACCTAATGGAGCATGCTATTCCTGTGTTAGTTATCACTTCGAAGGGGGTGCAAGAGTGTAAGTAATGGGT

'''
        testMaf = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test.maf')),
                               block, g_headers)
        parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        cmd = [os.path.abspath(os.path.join(parent, 'test', 'mafBlockExtractor'))]
        cmd += ['--maf', os.path.abspath(os.path.join(tmpDir, 'test.maf')),
                '--seq', 'simHuman.chrA', '--start', '%d' % 29953315, 
                '--stop', '%d' % 29953315,
                '--soft']
        outpipes = [os.path.abspath(os.path.join(tmpDir, 'extracted.maf'))]
        mtt.recordCommands([cmd], tmpDir, outPipes=outpipes)
        mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
        self.assertTrue(mtt.mafIsEmpty(os.path.join(tmpDir, 'extracted.maf'), g_headers))
        self.assertTrue(mafval.validateMaf(os.path.join(tmpDir, 'extracted.maf'), customOpts))
        mtt.removeDir(tmpDir)
    def testMemory0(self):
        """ If valgrind is installed on the system, check for memory related errors (0).
        """
        valgrind = mtt.which('valgrind')
        if valgrind is None:
            return
        mtt.makeTempDirParent()
        for i in xrange(0, 10):
            shuffledBlocks = []
            tmpDir = os.path.abspath(mtt.makeTempDir('memory_0'))
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
            testMaf = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test_%d.maf' % i)),
                                   ''.join(shuffledBlocks), g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = mtt.genericValgrind(tmpDir)
            cmd.append(os.path.abspath(os.path.join(parent, 'test', 'mafBlockExtractor')))
            cmd += ['--maf', os.path.abspath(os.path.join(tmpDir, 'test_%d.maf' % i)),
                    '--seq', g_targetSeq, '--start', '%d' % g_targetRange[0], 
                    '--stop', '%d' % g_targetRange[1],
                    '--soft']
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'extracted.maf'))]
            mtt.recordCommands([cmd], tmpDir, outPipes=outpipes)
            mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
            self.assertTrue(mtt.noMemoryErrors(os.path.join(tmpDir, 'valgrind.xml')))
            mtt.removeDir(tmpDir)
    def testMemory1(self):
        """ If valgrind is installed on the system, check for memory related errors (1).
        """
        valgrind = mtt.which('valgrind')
        if valgrind is None:
            return
        mtt.makeTempDirParent()
        for i in xrange(0, 10):
            tmpDir = os.path.abspath(mtt.makeTempDir('memory_1'))
            random.shuffle(g_nonOverlappingBlocks)
            testMaf = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test_%d.maf' % i)),
                                   ''.join(g_nonOverlappingBlocks), g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = mtt.genericValgrind(tmpDir)
            cmd.append(os.path.abspath(os.path.join(parent, 'test', 'mafBlockExtractor')))
            cmd += ['--maf', os.path.abspath(os.path.join(tmpDir, 'test_%d.maf' % i)),
                    '--seq', g_targetSeq, '--start', '%d' % g_targetRange[0], 
                    '--stop', '%d' % g_targetRange[1],
                    '--soft']
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'extracted.maf'))]
            mtt.recordCommands([cmd], tmpDir, outPipes=outpipes)
            mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
            self.assertTrue(mtt.noMemoryErrors(os.path.join(tmpDir, 'valgrind.xml')))
            mtt.removeDir(tmpDir)
class MafBlock:
    def __init__(self, aLine=None, mafLineList=None):
        self.aLine = aLine
        if mafLineList is None:
            self.mafLines = []
        else:
            self.mafLines = mafLineList
    def __str__(self):
        v = '%s\n' % self.aLine
        for ml in self.mafLines:
            v += '%s\n' % str(ml)
        v += '\n'
        return v
    def hashify(self):
        return self.__str__().replace(' ', '').replace('\n', '')
class MafLine:
    def __init__(self, name, start, length, strand, sourceLength, sequence):
        self.name = name
        self.start = start
        self.length = length
        self.strand = strand
        self.sourceLength = sourceLength
        self.sequence = sequence
    def __str__(self):
        return ('s %-15s %3d %3d %s %12d %s' 
                % (self.name, self.start, self.length, self.strand, self.sourceLength, self.sequence))
def rawBlockToObj(block):
    block = block.split('\n')
    assert(block[0].startswith('a'))
    b = MafBlock(block[0])
    for line in block[1:]:
        line = line.strip()
        if line == '':
            continue
        d = line.split()
        b.mafLines.append(MafLine(d[1], int(d[2]), int(d[3]), d[4], int(d[5]), d[6]))
    return b
class HardTrimTest(unittest.TestCase):
    testSet = [ # test0
        ([MafBlock('a score=0 ', [MafLine('target.chr0', 0, 13, '+', 158545518, 'gcagctgaaaaca'),
                                 MafLine('name.chr1', 0, 10, '+', 100, 'ATGT---ATGCCG'),
                                 MafLine('name2.chr1', 0, 10, '+', 100, 'ATGT---ATGCCG'),
                                 ])],
         'target.chr0',
         0,
         100,
         [MafBlock('a score=0', [MafLine('target.chr0', 0, 13, '+', 158545518, 'gcagctgaaaaca'),
                                 MafLine('name.chr1', 0, 10, '+', 100, 'ATGT---ATGCCG'),
                                 MafLine('name2.chr1', 0, 10, '+', 100, 'ATGT---ATGCCG'),
                                 ])],
         ),
        # test1
        ([MafBlock('a score=0', [MafLine('target.chr0', 0, 13, '+', 158545518, 'gcagctgaaaaca'),
                                 MafLine('name.chr1', 0, 10, '+', 100, 'ATGT---ATGCCG'),
                                 MafLine('name2.chr1', 0, 10, '+', 100, 'ATGT---ATGCCG'),
                                 ])],
         'target.chr0',
         2,
         20,
         [MafBlock('a score=0 mafBlockExtractor_splicedBlock=true splice_id=0_0',
                   [MafLine('target.chr0', 2, 11, '+', 158545518, 'agctgaaaaca'),
                    MafLine('name.chr1', 2, 8, '+', 100, 'GT---ATGCCG'),
                    MafLine('name2.chr1', 2, 8, '+', 100, 'GT---ATGCCG'),
                    ])],
         ),
        # test2
        ([MafBlock('a score=0', [MafLine('target.chr0', 0, 13, '+', 158545518, 'gcagctgaaaaca'),
                                 MafLine('name.chr1', 0, 10, '+', 100, 'ATGT---ATGCCG'),
                                 MafLine('name2.chr1', 0, 10, '+', 100, 'ATGT---ATGCCG'),
                                 ])],
         'target.chr0',
         0,
         10,
         [MafBlock('a score=0 mafBlockExtractor_splicedBlock=true splice_id=0_0',
                   [MafLine('target.chr0', 0, 11, '+', 158545518, 'gcagctgaaaa'),
                    MafLine('name.chr1', 0, 8, '+', 100, 'ATGT---ATGC'),
                    MafLine('name2.chr1', 0, 8, '+', 100, 'ATGT---ATGC'),
                    ])],
         ),
        # test3
        ([MafBlock('a score=0', [MafLine('target.chr0', 0, 13, '+', 158545518, 'gcagctgaaaaca'),
                                 MafLine('name.chr1', 0, 10, '+', 100, 'ATGT---ATGCCG'),
                                 MafLine('name2.chr1', 0, 10, '+', 100, 'ATGT---ATGCCG'),
                                 ])],
         'target.chr0',
         2,
         10,
         [MafBlock('a score=0 mafBlockExtractor_splicedBlock=true splice_id=0_0',
                   [MafLine('target.chr0', 2, 9, '+', 158545518, 'agctgaaaa'),
                    MafLine('name.chr1', 2, 6, '+', 100, 'GT---ATGC'),
                    MafLine('name2.chr1', 2, 6, '+', 100, 'GT---ATGC'),
                    ])],
         ),
        # test4
        ([MafBlock('a score=0', [MafLine('target.chr0', 0, 13, '+', 158545518, 'gcagctgaa---aaca'),
                                 MafLine('name.chr1', 0, 13, '+', 100,         'ATGT---ATGTAGCCG'),
                                 MafLine('name2.chr1', 0, 13, '+', 100,        'ATGT---ATGTAGCCG'),
                                 ])],
         'target.chr0',
         2,
         10,
         [MafBlock('a score=0 mafBlockExtractor_splicedBlock=true splice_id=0_0',
                   [MafLine('target.chr0', 2, 7, '+', 158545518, 'agctgaa'),
                    MafLine('name.chr1', 2, 4, '+', 100,         'GT---AT'),
                    MafLine('name2.chr1', 2, 4, '+', 100,        'GT---AT'),
                    ]),
          MafBlock('a score=0 mafBlockExtractor_splicedBlock=true splice_id=0_1', 
                   [MafLine('target.chr0', 9, 2, '+', 158545518, 'aa'),
                    MafLine('name.chr1', 9, 2, '+', 100,         'GC'),
                    MafLine('name2.chr1', 9, 2, '+', 100,        'GC'),
                    ]),
          ],
         ),
        # test5
        ([MafBlock('a score=0', [MafLine('target.chr0', 0, 13, '+', 158545518, 'gc-agc--tg-a-a---aaca'),
                                 MafLine('name.chr1', 0, 16, '+', 100,         'ATTGT-----AAGTGTAGCCG'),
                                 MafLine('name2.chr1', 0, 16, '+', 100,        'ATTGT-----AAGTGTAGCCG'),
                                 MafLine('name3.chr2', 0, 21, '+', 100,        'ATTGTAGTTTAAGTGTAGCCG'),
                                 ])],
         'target.chr0',
         2,
         10,
         [MafBlock('a score=0 mafBlockExtractor_splicedBlock=true splice_id=0_0',
                   [MafLine('target.chr0', 2, 3, '+', 158545518, 'agc'),
                    MafLine('name.chr1', 3, 2, '+', 100,         'GT-'),
                    MafLine('name2.chr1', 3, 2, '+', 100,        'GT-'),
                    MafLine('name3.chr2', 3, 3, '+', 100,        'GTA'),
                    ]),
          MafBlock('a score=0 mafBlockExtractor_splicedBlock=true splice_id=0_1',
                   [MafLine('target.chr0', 5, 2, '+', 158545518, 'tg'),
                    MafLine('name3.chr2', 8, 2, '+', 100,        'TT'),
                    ]),
          MafBlock('a score=0 mafBlockExtractor_splicedBlock=true splice_id=0_2', 
                   [MafLine('target.chr0', 7, 1, '+', 158545518, 'a'),
                    MafLine('name.chr1', 6, 1, '+', 100,         'A'),
                    MafLine('name2.chr1', 6, 1, '+', 100,        'A'),
                    MafLine('name3.chr2', 11, 1, '+', 100,       'A'),
                    ]),
          MafBlock('a score=0 mafBlockExtractor_splicedBlock=true splice_id=0_3', 
                   [MafLine('target.chr0', 8, 1, '+', 158545518, 'a'),
                    MafLine('name.chr1', 8, 1, '+', 100,         'T'),
                    MafLine('name2.chr1', 8, 1, '+', 100,        'T'),
                    MafLine('name3.chr2', 13, 1, '+', 100,       'T'),
                    ]),
          MafBlock('a score=0 mafBlockExtractor_splicedBlock=true splice_id=0_4', 
                   [MafLine('target.chr0', 9, 2, '+', 158545518, 'aa'),
                    MafLine('name.chr1', 12, 2, '+', 100,        'GC'),
                    MafLine('name2.chr1', 12, 2, '+', 100,       'GC'),
                    MafLine('name3.chr2', 17, 2, '+', 100,       'GC'),
                    ]),
          ],
         ),
        # test6
        ([MafBlock('a score=0', [MafLine('target.chr0', 0, 13, '-', 20, 'gcagctgaaaaca'),
                                 MafLine('name.chr1', 0, 10, '+', 100,  'ATTGT---AAGTG'),
                                 MafLine('name2.chr1', 0, 10, '+', 100, 'ATTGT---AAGTG'),
                                 MafLine('name3.chr2', 0, 9, '+', 100,  'ATTGT----AGTG'),
                                 ]),
          ],
         'target.chr0',
         0,
         19,
         [MafBlock('a score=0', [MafLine('target.chr0', 0, 13, '-', 20, 'gcagctgaaaaca'),
                                 MafLine('name.chr1', 0, 10, '+', 100,  'ATTGT---AAGTG'),
                                 MafLine('name2.chr1', 0, 10, '+', 100, 'ATTGT---AAGTG'),
                                 MafLine('name3.chr2', 0, 9, '+', 100,  'ATTGT----AGTG'),
                                 ]),
          ],
         ),
        # test7
        ([MafBlock('a score=0', [MafLine('target.chr0', 0, 13, '-', 20, 'gcagctgaaaaca'),
                                 MafLine('name.chr1', 0, 10, '+', 100,  'ATTGT---AAGTG'),
                                 MafLine('name2.chr1', 0, 10, '+', 100, 'ATTGT---AAGTG'),
                                 MafLine('name3.chr2', 0, 9, '+', 100,  'ATTGT----AGTG'),
                                 ]),
          ],
         'target.chr0',
         6,
         19,
         [MafBlock('a score=0', [MafLine('target.chr0', 0, 13, '-', 20, 'gcagctgaaaaca'),
                                 MafLine('name.chr1', 0, 10, '+', 100,  'ATTGT---AAGTG'),
                                 MafLine('name2.chr1', 0, 10, '+', 100, 'ATTGT---AAGTG'),
                                 MafLine('name3.chr2', 0, 9, '+', 100,  'ATTGT----AGTG'),
                                 ]),
          ],
         ),
        ####################
        # test8
        ([MafBlock('a score=0', [MafLine('target.chr0', 0, 13, '-',  20, 'gcagctgaaaaca'),
                                 MafLine('name.chr1',   0, 10, '+', 100, 'ATTGT---AAGTG'),
                                 MafLine('name2.chr1',  0, 10, '+', 100, 'ATTGT---AAGTG'),
                                 MafLine('name3.chr2',  0,  9, '+', 100, 'ATTGT----AGTG'),
                                 ]),
          ],
         'target.chr0',
         8,
         16,
         [MafBlock('a score=0 mafBlockExtractor_splicedBlock=true  splice_id=0_0',
                   [MafLine('target.chr0', 3, 9, '-',  20, 'gctgaaaac'),
                    MafLine('name.chr1',   3, 6, '+', 100, 'GT---AAGT'),
                    MafLine('name2.chr1',  3, 6, '+', 100, 'GT---AAGT'),
                    MafLine('name3.chr2',  3, 5, '+', 100, 'GT----AGT'),
                    ]),
          ],
         ),
        ####################
        # test9
        ([MafBlock('a score=0', [MafLine('target.chr0', 0, 13, '-', 20, 'g-c-a-g-c-t-g-a-a-a-a-c-a-'),
                                 MafLine('name.chr1', 0, 23, '+', 100,  'ATTGT---AAGTGTTTTTTTTTTTTT'),
                                 MafLine('name2.chr1', 0, 23, '+', 100, 'ATTGT---AAGTGTTTTTTTTTTTTT'),
                                 MafLine('name3.chr2', 0, 22, '+', 100, 'ATTGT----AGTGTTTTTTTTTTTTT'),
                                 ]),
          ],
         'target.chr0',
         0,
         20,
         [MafBlock('a score=0 mafBlockExtractor_splicedBlock=true splice_id=0_0', 
                   [MafLine('target.chr0', 0, 1, '-', 20, 'g'),
                    MafLine('name.chr1',  0, 1, '+', 100, 'A'),
                    MafLine('name2.chr1', 0, 1, '+', 100, 'A'),
                    MafLine('name3.chr2', 0, 1, '+', 100, 'A'),
                    ]),
          MafBlock('a score=0 mafBlockExtractor_splicedBlock=true  splice_id=0_1', 
                   [MafLine('target.chr0', 1, 1, '-', 20, 'c'),
                    MafLine('name.chr1',  2, 1, '+', 100, 'T'),
                    MafLine('name2.chr1', 2, 1, '+', 100, 'T'),
                    MafLine('name3.chr2', 2, 1, '+', 100, 'T'),
                    ]),
          MafBlock('a score=0 mafBlockExtractor_splicedBlock=true splice_id=0_2', 
                   [MafLine('target.chr0', 2, 1, '-', 20, 'a'),
                    MafLine('name.chr1',  4, 1, '+', 100, 'T'),
                    MafLine('name2.chr1', 4, 1, '+', 100, 'T'),
                    MafLine('name3.chr2', 4, 1, '+', 100, 'T'),
                    ]),
          MafBlock('a score=0 mafBlockExtractor_splicedBlock=true splice_id=0_3', 
                   [MafLine('target.chr0', 3, 1, '-', 20, 'g'),
                    ]),
          MafBlock('a score=0 mafBlockExtractor_splicedBlock=true splice_id=0_4', 
                   [MafLine('target.chr0', 4, 1, '-', 20, 'c'),
                    MafLine('name.chr1',  5, 1, '+', 100, 'A'),
                    MafLine('name2.chr1', 5, 1, '+', 100, 'A'),
                    ]),
          MafBlock('a score=0 mafBlockExtractor_splicedBlock=true splice_id=0_5', 
                   [MafLine('target.chr0', 5, 1, '-', 20, 't'),
                    MafLine('name.chr1',  7, 1, '+', 100, 'G'),
                    MafLine('name2.chr1', 7, 1, '+', 100, 'G'),
                    MafLine('name3.chr2', 6, 1, '+', 100, 'G'),
                    ]),
          MafBlock('a score=0 mafBlockExtractor_splicedBlock=true splice_id=0_6', 
                   [MafLine('target.chr0', 6, 1, '-',  20, 'g'),
                    MafLine('name.chr1',   9, 1, '+', 100, 'G'),
                    MafLine('name2.chr1',  9, 1, '+', 100, 'G'),
                    MafLine('name3.chr2',  8, 1, '+', 100, 'G'),
                    ]),
          MafBlock('a score=0 mafBlockExtractor_splicedBlock=true splice_id=0_7', 
                   [MafLine('target.chr0', 7, 1, '-',  20, 'a'),
                    MafLine('name.chr1',  11, 1, '+', 100, 'T'),
                    MafLine('name2.chr1', 11, 1, '+', 100, 'T'),
                    MafLine('name3.chr2', 10, 1, '+', 100, 'T'),
                    ]),
          MafBlock('a score=0 mafBlockExtractor_splicedBlock=true splice_id=0_8', 
                   [MafLine('target.chr0', 8, 1, '-',  20, 'a'),
                    MafLine('name.chr1',  13, 1, '+', 100, 'T'),
                    MafLine('name2.chr1', 13, 1, '+', 100, 'T'),
                    MafLine('name3.chr2', 12, 1, '+', 100, 'T'),
                    ]),
          MafBlock('a score=0 mafBlockExtractor_splicedBlock=true splice_id=0_9', 
                   [MafLine('target.chr0', 9, 1, '-',  20, 'a'),
                    MafLine('name.chr1',  15, 1, '+', 100, 'T'),
                    MafLine('name2.chr1', 15, 1, '+', 100, 'T'),
                    MafLine('name3.chr2', 14, 1, '+', 100, 'T'),
                    ]),
          MafBlock('a score=0 mafBlockExtractor_splicedBlock=true splice_id=0_10', 
                   [MafLine('target.chr0', 10, 1, '-',  20, 'a'),
                    MafLine('name.chr1',   17, 1, '+', 100, 'T'),
                    MafLine('name2.chr1',  17, 1, '+', 100, 'T'),
                    MafLine('name3.chr2',  16, 1, '+', 100, 'T'),
                    ]),
          MafBlock('a score=0 mafBlockExtractor_splicedBlock=true splice_id=0_11', 
                   [MafLine('target.chr0', 11, 1, '-',  20, 'c'),
                    MafLine('name.chr1',   19, 1, '+', 100, 'T'),
                    MafLine('name2.chr1',  19, 1, '+', 100, 'T'),
                    MafLine('name3.chr2',  18, 1, '+', 100, 'T'),
                    ]),
          MafBlock('a score=0 mafBlockExtractor_splicedBlock=true splice_id=0_12', 
                   [MafLine('target.chr0', 12, 1, '-',  20, 'a'),
                    MafLine('name.chr1',   21, 1, '+', 100, 'T'),
                    MafLine('name2.chr1',  21, 1, '+', 100, 'T'),
                    MafLine('name3.chr2',  20, 1, '+', 100, 'T'),
                    ]),
          ],
         ),
        # test10
        ([MafBlock('a score=0', [MafLine('target.chr0', 0, 13, '+', 158545518, 'gcagctgaaaaca'),
                                 MafLine('name.chr1', 0, 10, '+', 100, 'ATGT---ATGCCG'),
                                 MafLine('name2.chr1', 0, 10, '+', 100, 'ATGT---ATGCCG'),
                                 ])],
         'target.chr0',
         20,
         30,
         [
          ],
         ),
        ]
    def testHardExtraction_0(self):
        """ mafBlockExtractor should output trimmed blocks by default.
        """
        valgrind = mtt.which('valgrind')
        if valgrind is None:
            return
        mtt.makeTempDirParent()
        customOpts = GenericValidationOptions()
        i = -1
        for inblocks, seq, start, stop, outblocks in self.testSet:
            i += 1
            tmpDir = os.path.abspath(mtt.makeTempDir('hardextraction_0'))
            testMaf = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test_%d.maf' % i)),
                                   ''.join(map(str, inblocks)), g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = [os.path.abspath(os.path.join(parent, 'test', 'mafBlockExtractor'))]
            cmd += ['--maf', os.path.abspath(os.path.join(tmpDir, 'test_%d.maf' % i)),
                    '--seq', seq, '--start', '%d' % start, 
                    '--stop', '%d' % stop,
                    ]
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'extracted.maf'))]
            mtt.recordCommands([cmd], tmpDir, outPipes=outpipes)
            mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
            self.assertTrue(mafIsHardExtracted('hardextraction_0_%d' % i, outblocks, 
                                               os.path.join(tmpDir, 'extracted.maf')))
            self.assertTrue(mafval.validateMaf(os.path.join(tmpDir, 'extracted.maf'), customOpts))
            mtt.removeDir(tmpDir)
    def testMemory2(self):
        """ If valgrind is installed on the system, check for memory related errors (2).
        """
        valgrind = mtt.which('valgrind')
        if valgrind is None:
            return
        mtt.makeTempDirParent()
        i = -1
        for inblocks, seq, start, stop, outblocks in self.testSet:
            i += 1
            tmpDir = os.path.abspath(mtt.makeTempDir('memory_2'))
            random.shuffle(g_nonOverlappingBlocks)
            testMaf = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test_%d.maf' % i)),
                                   ''.join(map(str, inblocks)), g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = mtt.genericValgrind(tmpDir)
            cmd.append(os.path.abspath(os.path.join(parent, 'test', 'mafBlockExtractor')))
            cmd += ['--maf', os.path.abspath(os.path.join(tmpDir, 'test_%d.maf' % i)),
                    '--seq', seq, '--start', '%d' % start, 
                    '--stop', '%d' % stop,
                    ]
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'extracted.maf'))]
            mtt.recordCommands([cmd], tmpDir, outPipes=outpipes)
            mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
            self.assertTrue(mtt.noMemoryErrors(os.path.join(tmpDir, 'valgrind.xml')))
            mtt.removeDir(tmpDir)

class CuTestMemory(unittest.TestCase):
    def test_CuTestMemory(self):
        """ If valgrind is installed on the system, check for memory related errors in CuTests.
        """
        mtt.makeTempDirParent()
        valgrind = mtt.which('valgrind')
        if valgrind is None:
            return
        tmpDir = os.path.abspath(mtt.makeTempDir('memory_CuTest'))
        parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        cmd = mtt.genericValgrind(tmpDir)
        cmd.append(os.path.abspath(os.path.join(parent, 'test', 'allTests')))
        outpipes = [os.path.join('/dev', 'null')]
        mtt.recordCommands([cmd], tmpDir)
        mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
        self.assertTrue(mtt.noMemoryErrors(os.path.join(tmpDir, 'valgrind.xml')))
        mtt.removeDir(tmpDir)

if __name__ == '__main__':
    unittest.main()
