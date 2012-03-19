import os
import random
import shutil
import subprocess
import sys
import unittest

target = 'hg18.chr7'

header = '''##maf version=1 scoring=tba.v8
# tba.v8 (((human chimp) baboon) (mouse rat))

'''

# does not contain target
nonTargetBlocks = ['''a score=23261.0
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon.chr0    249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA

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

# target is hg18.chr7, blocks are in correct order.
targetBlocks = ['''a score=23263.0
s hg18.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon.chr0    116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
                   
''',
                '''a score=5062.0                    
s hg18.chr7    27699739 6 + 158545518 TAAAGA
s panTro1.chr6 28862317 6 + 161576975 TAAAGA
s baboon.chr0    241163 6 +   4622798 TAAAGA 
s mm4.chr6     53303881 6 + 151104725 TAAAGA
s rn3.chr4     81444246 6 + 187371129 taagga

''',
                '''a score=23264.0
s hg18.chr7    27707000 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon.chr0    249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA

''',
                '''a score=6636.0
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon.chr0    249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA
s hg18.chr7    27707221 13 + 158545518 gcagctgaaaaca

''',
                '''a score=0
s hg18.chr7              237249719 26 - 247249719 TTTTTGAAAAACAAACAACAAGTTGG
s panTro2.chrUn            9697231 26 +  58616431 TTTTTGAAAAACAAACAACAAGTTGG
q panTro2.chrUn                                   99999999999999999999999999
s dasNov1.scaffold_179265     1474  7 +      4584 TT----------AAGCA---------
q dasNov1.scaffold_179265                         99----------32239--------- 

''',
                ]
def testFile(s):
    makeTempDir()
    mafFile = os.path.abspath(os.path.join(os.curdir, 'tempTestDir', 'test.maf'))
    f = open(mafFile, 'w')
    f.write(header)
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
def mafIsSorted(maf):
    f = open(maf)
    # header lines
    l = f.next()
    l = f.next()
    l = f.next()
    nonTBlocks = []
    expectedNonTargetBlocks = list(nonTargetBlocks)
    expectedTargetBlocks = list(targetBlocks)
    for i in xrange(0, len(expectedNonTargetBlocks)):
        nonTBlocks.append(extractBlockStr(f).replace(' ', ''))
        expectedNonTargetBlocks[i] = expectedNonTargetBlocks[i].replace(' ', '')
    if nonTBlocks != expectedNonTargetBlocks:
        print '\nexpected'
        print expectedNonTargetBlocks
        print 'observed'
        print nonTBlocks

        f.close()
        return False
    sortedBlocks = []
    for i in xrange(0, len(expectedTargetBlocks)):
        sortedBlocks.append(extractBlockStr(f).replace(' ', ''))
        expectedTargetBlocks[i] = expectedTargetBlocks[i].replace(' ', '')
    if sortedBlocks != expectedTargetBlocks:
        f.close()
        return False
    f.close()
    return True
def extractBlockStr(f):
    block = ''
    for line in f:
        line = line.strip()
        if line == '':
            return block + '\n'
        block += '%s\n' % line
    return block + '\n'
    
class SortTest(unittest.TestCase):
    def testSorting(self):
        """ Blocks should be sorted by the start field of the target sequence,
        blocks that do not contain the target sequence should appear in the 
        output at the start of the file, in the same order they appear in the input.
        """
        shuffledTargets = list(targetBlocks)
        for i in xrange(0, 100):
            tmpDir = os.path.abspath(makeTempDir())
            random.shuffle(nonTargetBlocks)
            random.shuffle(shuffledTargets)
            shuffledBlocks = list(shuffledTargets)
            lower = 0
            for j in xrange(0, len(nonTargetBlocks)):
                # randomly insert the non target blocks, but keep a record
                # of their relative order.
                index = random.randint(lower, len(shuffledBlocks))
                shuffledBlocks.insert(index, nonTargetBlocks[j])
                lower = index + 1
            testMaf = testFile(''.join(shuffledBlocks))
            binParent = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
            cmd = [os.path.abspath(os.path.join(binParent, 'bin', 'mafBlockSorter'))]
            cmd += ['--seq', 'hg18.chr7']
            inpipes = [testMaf]
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'sorted.maf'))]
            runCommandsS([cmd], tmpDir, inPipes=inpipes, outPipes=outpipes)
            self.assertTrue(mafIsSorted(os.path.join(tmpDir, 'sorted.maf')))
            removeTempDir()

if __name__ == '__main__':
    unittest.main()
