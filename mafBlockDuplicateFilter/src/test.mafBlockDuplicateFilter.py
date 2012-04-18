import os
import random
import shutil
import subprocess
import sys
import unittest

g_header = '''##maf version=1 scoring=tba.v8
# tba.v8 (((human chimp) baboon) (mouse rat))

'''

g_duplicateBlocks = [('''a score=0
#dup block 1
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
#dup block 1
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
#dup block 2
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
#dup block 2
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
#dup block 3
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
#dup block 3
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
#dup block 4
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
#dup block 4
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

def testFile(s):
    makeTempDir()
    mafFile = os.path.abspath(os.path.join(os.curdir, 'tempTestDir', 'test.maf'))
    f = open(mafFile, 'w')
    f.write(g_header)
    f.write(s)
    f.close()
    return mafFile
def makeTempDir():
    if not os.path.exists(os.path.join(os.curdir, 'tempTestDir')):
        os.mkdir(os.path.join(os.curdir, 'tempTestDir'))
    return os.path.join(os.curdir, 'tempTestDir')
def removeTempDir():
    if os.path.exists(os.path.join(os.curdir, 'tempTestDir')):
        shutil.rmtree(os.path.join(os.curdir, 'tempTestDir'))
def runCommandsS(cmds, localTempDir, inPipes=[], outPipes=[]):
    """ 
    runCommandsS uses the subprocess module
    to issue serial processes from the cmds list.
    """
    if not len(inPipes):
        inPipes = [None] * len(cmds)
    if not len(outPipes):
        outPipes = [None] * len(cmds)
    for i, c in enumerate(cmds, 0):
        if inPipes[i] is None:
            sin = None
        else:
            sin = subprocess.PIPE
        if outPipes[i] is None:
            sout = None
        else:
            sout = subprocess.PIPE
        p = subprocess.Popen(c, cwd=localTempDir, stdin=sin, stdout=sout)
            
        if inPipes[i] is None:
            sin = None
        else:
            if not os.path.exists(inPipes[i]):
                raise IOError('Unable to locate inPipe file: %s for command %s' % (inPipes[i], ' '.join(c)))
            sin = open(inPipes[i], 'r').read()
        if outPipes[i] is None:
            pout, perr = p.communicate(sin)
            handleReturnCode(p.returncode, cmds[i])
        else:
            f = open(outPipes[i], 'w')
            f.write(p.communicate(sin)[0])
            f.close()
            handleReturnCode(p.returncode, cmds[i])
def handleReturnCode(retcode, cmd):
    if not isinstance(retcode, int):
        raise TypeError('handleReturnCode takes an integer for '
                        'retcode, not a %s.' % retcode.__class__)
    if not isinstance(cmd, list):
        raise TypeError('handleReturnCode takes a list for '
                        'cmd, not a %s.' % cmd.__class__)
    if retcode:
        if retcode < 0:
            raise RuntimeError('Experienced an error while trying to execute: '
                               '%s SIGNAL:%d' %(' '.join(cmd), -retcode))
        else:
            raise RuntimeError('Experienced an error while trying to execute: '
                               '%s retcode:%d' %(' '.join(cmd), retcode))
def mafIsEmpty(maf):
    f = open(maf)
    s = f.read()
    if s == g_header:
        return True
    return False
def mafIsFiltered(maf, blockList):
    f = open(maf)
    # header lines
    l = f.next()
    l = f.next()
    l = f.next()
    for i in xrange(0, len(blockList)):
        # walk through the maf, assessing the equavalence to the blockList items
        b = extractBlockStr(f)
        if b != blockList[i]:
            print 'extracted block '
            print b
            print 'expected block '
            print blockList[i]
            return False
    return True
def extractBlockStr(f):
    block = ''
    for line in f:
        line = line.strip()
        if line == '':
            return block + '\n'
        block += '%s\n' % line
    return block + '\n'
    
class ExtractionTest(unittest.TestCase):
    def testFilter(self):
        """ mafBlockDuplicateFilter should filter out duplicates in blocks according to sequence similarity to the consensus.
        """
        for i in xrange(0, 10):
            shuffledBlocks = []
            expectedOutput = []
            tmpDir = os.path.abspath(makeTempDir())
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
            testMaf = testFile(''.join(shuffledBlocks))
            binParent = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
            cmd = [os.path.abspath(os.path.join(binParent, 'bin', 'mafBlockDuplicateFilter'))]
            inpipes = [testMaf]
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'filtered.maf'))]
            runCommandsS([cmd], tmpDir, inPipes=inpipes, outPipes=outpipes)
            self.assertTrue(mafIsFiltered(os.path.join(tmpDir, 'filtered.maf'), expectedOutput))
            removeTempDir()
    def testNonExtraction(self):
        """ mafBlockExtractor should not filter out any sequences from blocks when there are no duplicates.
        """
        for i in xrange(0, 10):
            tmpDir = os.path.abspath(makeTempDir())
            random.shuffle(g_nonDuplicateBlocks)
            testMaf = testFile(''.join(g_nonDuplicateBlocks))
            expectedOutput = g_nonDuplicateBlocks
            binParent = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
            cmd = [os.path.abspath(os.path.join(binParent, 'bin', 'mafBlockDuplicateFilter'))]
            inpipes = [testMaf]
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'filtered.maf'))]
            runCommandsS([cmd], tmpDir, inPipes=inpipes, outPipes=outpipes)
            self.assertTrue(mafIsFiltered(os.path.join(tmpDir, 'filtered.maf'), expectedOutput))
            removeTempDir()

if __name__ == '__main__':
    unittest.main()
