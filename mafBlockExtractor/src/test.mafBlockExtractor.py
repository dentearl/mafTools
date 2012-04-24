import os
import random
import shutil
import subprocess
import sys
import unittest
import xml.etree.ElementTree as ET
import xml.parsers.expat

g_targetSeq = 'target.chr0'
g_targetRange = (50, 70) # zero based, inclusive

g_header = '''##maf version=1 scoring=tba.v8
# tba.v8 (((human chimp) baboon) (mouse rat))

'''

g_overlappingBlocks = ['''a score=0
s target.chr0        38 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA
s name.chr1           0 10 +       100 ATGT---ATGCCG
s name2.chr1         50 10 +       100 ATGT---ATGCCG
s name3.chr9         50 10 +       100 ATGTA---TGCCG
s name4.chr&         50 10 +       100 ATG---TATGCCG
s name5 50 10 + 100 ATGTATGCCG

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
def mafIsExtracted(maf):
    f = open(maf)
    # three header lines
    l = f.next()
    l = f.next()
    l = f.next()
    for i in xrange(0, len(g_overlappingBlocks)):
        b = extractBlockStr(f)
        if b not in g_overlappingBlocks:
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
    def testExtraction(self):
        """ mafBlockExtractor should output blocks that meet the criteria for extraction. That is they contain the taget sequence and have at least one base in the target range.
        """
        for i in xrange(0, 10):
            shuffledBlocks = []
            tmpDir = os.path.abspath(makeTempDir())
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
            testMaf = testFile(''.join(shuffledBlocks))
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = [os.path.abspath(os.path.join(parent, 'test', 'mafBlockExtractor'))]
            cmd += ['--seq', g_targetSeq, '--start', '%d' % g_targetRange[0], 
                    '--stop', '%d' % g_targetRange[1]]
            inpipes = [testMaf]
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'extracted.maf'))]
            runCommandsS([cmd], tmpDir, inPipes=inpipes, outPipes=outpipes)
            self.assertTrue(mafIsExtracted(os.path.join(tmpDir, 'extracted.maf')))
            removeTempDir()
    def testNonExtraction(self):
        """ mafBlockExtractor should not extract blocks when they do not match.
        """
        for i in xrange(0, 10):
            tmpDir = os.path.abspath(makeTempDir())
            random.shuffle(g_nonOverlappingBlocks)
            testMaf = testFile(''.join(g_nonOverlappingBlocks))
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = [os.path.abspath(os.path.join(parent, 'test', 'mafBlockExtractor'))]
            cmd += ['--seq', g_targetSeq, '--start', '%d' % g_targetRange[0], 
                    '--stop', '%d' % g_targetRange[1]]
            inpipes = [testMaf]
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'extracted.maf'))]
            runCommandsS([cmd], tmpDir, inPipes=inpipes, outPipes=outpipes)
            self.assertTrue(mafIsEmpty(os.path.join(tmpDir, 'extracted.maf')))
            removeTempDir()
    def testMemory1(self):
        """ If valgrind is installed on the system, check for memory related errors (1).
        """
        valgrind = which('valgrind')
        if valgrind is None:
            return
        for i in xrange(0, 10):
            shuffledBlocks = []
            tmpDir = os.path.abspath(makeTempDir())
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
            testMaf = testFile(''.join(shuffledBlocks))
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = [valgrind, '--leak-check=yes', '--track-origins=yes', '--xml=yes', 
                   '--xml-file=' + os.path.join(tmpDir, 'valgrind.xml')]
            cmd.append(os.path.abspath(os.path.join(parent, 'test', 'mafBlockExtractor')))
            cmd += ['--seq', g_targetSeq, '--start', '%d' % g_targetRange[0], 
                    '--stop', '%d' % g_targetRange[1]]
            inpipes = [testMaf]
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'extracted.maf'))]
            runCommandsS([cmd], tmpDir, inPipes=inpipes, outPipes=outpipes)
            self.assertTrue(noMemoryErrors(os.path.join(tmpDir, 'valgrind.xml')))
            removeTempDir()
    def testMemory2(self):
        """ If valgrind is installed on the system, check for memory related errors (2).
        """
        valgrind = which('valgrind')
        if valgrind is None:
            return
        for i in xrange(0, 10):
            tmpDir = os.path.abspath(makeTempDir())
            random.shuffle(g_nonOverlappingBlocks)
            testMaf = testFile(''.join(g_nonOverlappingBlocks))
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = [valgrind, '--leak-check=yes', '--track-origins=yes', '--xml=yes', 
                   '--xml-file=' + os.path.join(tmpDir, 'valgrind.xml')]
            cmd.append(os.path.abspath(os.path.join(parent, 'test', 'mafBlockExtractor')))
            cmd += ['--seq', g_targetSeq, '--start', '%d' % g_targetRange[0], 
                    '--stop', '%d' % g_targetRange[1]]
            inpipes = [testMaf]
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'extracted.maf'))]
            runCommandsS([cmd], tmpDir, inPipes=inpipes, outPipes=outpipes)
            self.assertTrue(noMemoryErrors(os.path.join(tmpDir, 'valgrind.xml')))
            removeTempDir()

if __name__ == '__main__':
    unittest.main()
