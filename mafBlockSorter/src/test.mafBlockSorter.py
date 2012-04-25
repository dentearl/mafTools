import os
import random
import shutil
import subprocess
import sys
import unittest
import xml.etree.ElementTree as ET
import xml.parsers.expat

g_target = 'hg18.chr7'
g_headers = ['''##maf version=1 scoring=tba.v8
# tba.v8 (((human chimp) baboon) (mouse rat))

''',
           '''##maf version=1 scoring=tba.v8
# tba.v8 (((human chimp) baboon) (mouse rat))
''',]

# does not contain target
g_nonTargetBlocks = ['''a score=23261.0
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
g_targetBlocks = ['''a score=23263.0
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
def testFile(s):
    makeTempDir()
    mafFile = os.path.abspath(os.path.join(os.curdir, 'tempTestDir', 'test.maf'))
    f = open(mafFile, 'w')
    # choose one of the headers at random
    f.write(random.choice(g_headers))
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
def which(program):
    """which() acts like the unix utility which, but is portable between os.
    If the program does not exist in the PATH then 'None' is returned. 
    """
    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath != '':
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None
def noMemoryErrors(xml):
    try:
        tree = ET.parse(xml)
    except xml.parsers.expat.ExpatError:
        raise RuntimeError('Input xml, %s is not a well formed xml document.' % xml)
    root = tree.getroot()
    errors = root.findall('error')
    if len(errors):
        return False
    errorcounts = root.find('errorcounts')
    if len(errorcounts):
        return False
    return True
def processHeader(f):
    for line in f:
        line = line.strip()
        if line == '':
            return None
        if line.startswith('a'):
            return line
def mafIsSorted(maf):
    f = open(maf)
    lastLine = processHeader(f)
    observedNonTargetBlocks = []
    expectedNonTargetBlocks = list(g_nonTargetBlocks)
    for i in xrange(0, len(expectedNonTargetBlocks)):
        observedNonTargetBlocks.append(extractBlockStr(f, lastLine).replace(' ', ''))
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
        sortedBlocks.append(extractBlockStr(f).replace(' ', ''))
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
def extractBlockStr(f, lastLine=None):
    if lastLine is None:
        block = ''
    else:
        block = lastLine + '\n'
    for line in f:
        line = line.strip()
        if line == '':
            return block + '\n'
        block += '%s\n' % line
    return block + '\n'
    
class SortTest(unittest.TestCase):
    def testSorting(self):
        """ Blocks should be sorted by the start field of the target sequence, blocks that do not contain the target sequence should appear in the output at the start of the file, in the same order they appear in the input.
        """
        shuffledTargets = list(g_targetBlocks)
        for i in xrange(0, 100):
            tmpDir = os.path.abspath(makeTempDir())
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
            testMaf = testFile(''.join(shuffledBlocks))
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = [os.path.abspath(os.path.join(parent, 'test', 'mafBlockSorter'))]
            cmd += ['--seq', 'hg18.chr7']
            inpipes = [testMaf]
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'sorted.maf'))]
            runCommandsS([cmd], tmpDir, inPipes=inpipes, outPipes=outpipes)
            self.assertTrue(mafIsSorted(os.path.join(tmpDir, 'sorted.maf')))
            removeTempDir()
    def testMemory1(self):
        """ If valgrind is installed on the system, check for memory related errors (1).
        """
        valgrind = which('valgrind')
        if valgrind is None:
            return
        shuffledTargets = list(g_targetBlocks)
        for i in xrange(0, 20):
            tmpDir = os.path.abspath(makeTempDir())
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
            testMaf = testFile(''.join(shuffledBlocks))
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = [valgrind, '--leak-check=yes', '--track-origins=yes', '--xml=yes', 
                   '--xml-file=' + os.path.join(tmpDir, 'valgrind.xml')]
            cmd.append(os.path.abspath(os.path.join(parent, 'test', 'mafBlockSorter')))
            cmd += ['--seq', 'hg18.chr7']
            inpipes = [testMaf]
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'sorted.maf'))]
            runCommandsS([cmd], tmpDir, inPipes=inpipes, outPipes=outpipes)
            self.assertTrue(noMemoryErrors(os.path.join(tmpDir, 'valgrind.xml')))
            removeTempDir()

if __name__ == '__main__':
    unittest.main()
