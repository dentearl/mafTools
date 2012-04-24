import os
import random
import re
import shutil
import subprocess
import sys
import unittest
import xml.etree.ElementTree as ET
import xml.parsers.expat

g_targetSeq = 'target.chr0'
g_targetRange = (50, 70) # zero based, inclusive
g_header = None
g_headers = ['''##maf version=1 scoring=tba.v8
# tba.v8 (((human chimp) baboon) (mouse rat))

''', 
             '''##maf version=1 scoring=tba.v8
# tba.v8 (((human chimp) baboon) (mouse rat))
''',]

# these are triples of blocks, target positions for True, 
# line numbers where the target occurs.
g_overlappingBlocks = [('''a score=0
s target.chr0        38 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA
s name.chr1           0 10 +       100 ATGT---ATGCCG
s name2.chr1         50 10 +       100 ATGT---ATGCCG
s name3.chr9         50 10 +       100 ATGTA---TGCCG
s name4.chr&         50 10 +       100 ATG---TATGCCG
s name5 50 10 + 100 ATGTATGCCG

''', 38, [2]),
            # the overlap is one base long, right on 50.
            ('''a score=0
s name                0 10 +       100 ATGTATGC---CG
s name2.chr1         50 10 +       100 ATGTATG---CCG
s name3.chr9         50 10 +       100 ATGTATG---CCG
s name4.chr&         50 10 +       100 ATGTATG---CCG
s name5              50 10 +       100 ATGTATG---CCG
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s target.chr0 158545457 10 - 158545518 ATGTATG---CCG

''', 60, [9]),
            # the overlap is 10 bases long, 51-61
            ('''a score=0
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

''', 70, [8]),
            # the overlap is 9 bases long, 62-70
            ]
g_nonOverlappingBlocks = [('''a score=23262.b0     
s hg18.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon         116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG

''', 0, None),
                          ('''a score=5062.0                    
s hg18.chr7    27699739 6 + 158545518 TAAAGA
s panTro1.chr6 28862317 6 + 161576975 TAAAGA
s baboon         241163 6 +   4622798 TAAAGA 
s mm4.chr6     53303881 6 + 151104725 TAAAGA
s rn3.chr4     81444246 6 + 187371129 taagga

''', 0, None),
                          ('''a score=6636.0
s hg18.chr7    27707221 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA

''', 0, None),
    ]

def testFile(s):
    global g_header
    makeTempDir()
    mafFile = os.path.abspath(os.path.join(os.curdir, 'tempTestDir', 'test.maf'))
    f = open(mafFile, 'w')
    g_header = random.choice(g_headers)
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
def foundLines(lineList, text):
    global g_header
    f = open(text, 'r')
    n = len(re.findall(re.compile('\n'), g_header))
    for line in f:
        line = line.strip()
        data = line.split(':')
        if int(data[0]) - n not in lineList:
            return False
    return True
def fileIsEmpty(filename):
    f = open(filename)
    s = f.read()
    if s == '':
        return True
    return False
    
class FindTest(unittest.TestCase):
    def testFind(self):
        """ mafBlockFinder should report information about matching sequences within blocks.
        """
        for i in xrange(0, len(g_overlappingBlocks)):
            tmpDir = os.path.abspath(makeTempDir())
            testMaf = testFile(g_overlappingBlocks[i][0])
            binParent = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
            cmd = [os.path.abspath(os.path.join(binParent, 'bin', 'mafBlockFinder'))]
            cmd += ['--seq', g_targetSeq, '--pos', '%d' % g_overlappingBlocks[i][1]]
            inpipes = [testMaf]
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'found.txt'))]
            runCommandsS([cmd], tmpDir, inPipes=inpipes, outPipes=outpipes)
            self.assertTrue(foundLines(g_overlappingBlocks[i][2], os.path.join(tmpDir, 'found.txt')))
            removeTempDir()
    def testNonFind(self):
        """ mafBlockFinder should not report any lines when blocks do not match.
        """
        for i in xrange(0, len(g_nonOverlappingBlocks)):
            tmpDir = os.path.abspath(makeTempDir())
            testMaf = testFile(''.join(g_nonOverlappingBlocks[i][0]))
            binParent = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
            cmd = [os.path.abspath(os.path.join(binParent, 'bin', 'mafBlockFinder'))]
            cmd += ['--seq', g_targetSeq, '--pos', '%d' % g_nonOverlappingBlocks[i][1]]
            inpipes = [testMaf]
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'found.txt'))]
            runCommandsS([cmd], tmpDir, inPipes=inpipes, outPipes=outpipes)
            self.assertTrue(fileIsEmpty(os.path.join(tmpDir, 'found.txt')))
            removeTempDir()
    def testMemory1(self):
        """ If valgrind is installed on the system, check for memory related errors (1).
        """
        valgrind = which('valgrind')
        if valgrind is None:
            return
        for i in xrange(0, len(g_overlappingBlocks)):
            tmpDir = os.path.abspath(makeTempDir())
            testMaf = testFile(g_overlappingBlocks[i][0])
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = [valgrind, '--leak-check=yes', '--track-origins=yes', '--xml=yes', 
                   '--xml-file=' + os.path.join(tmpDir, 'valgrind.xml')]
            cmd.append(os.path.abspath(os.path.join(parent, 'test', 'mafBlockFinder')))
            cmd += ['--seq', g_targetSeq, '--pos', '%d' % g_overlappingBlocks[i][1]]
            inpipes = [testMaf]
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'found.txt'))]
            runCommandsS([cmd], tmpDir, inPipes=inpipes, outPipes=outpipes)
            self.assertTrue(noMemoryErrors(os.path.join(tmpDir, 'valgrind.xml')))
            removeTempDir()
    def testMemory2(self):
        """ If valgrind is installed on the system, check for memory related errors (2).
        """
        valgrind = which('valgrind')
        if valgrind is None:
            return
        for i in xrange(0, len(g_nonOverlappingBlocks)):
            tmpDir = os.path.abspath(makeTempDir())
            testMaf = testFile(''.join(g_nonOverlappingBlocks[i][0]))
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = [valgrind, '--leak-check=yes', '--track-origins=yes', '--xml=yes', 
                   '--xml-file=' + os.path.join(tmpDir, 'valgrind.xml')]
            cmd.append(os.path.abspath(os.path.join(parent, 'test', 'mafBlockFinder')))
            cmd += ['--seq', g_targetSeq, '--pos', '%d' % g_nonOverlappingBlocks[i][1]]
            inpipes = [testMaf]
            outpipes = [os.path.abspath(os.path.join(tmpDir, 'found.txt'))]
            runCommandsS([cmd], tmpDir, inPipes=inpipes, outPipes=outpipes)
            self.assertTrue(noMemoryErrors(os.path.join(tmpDir, 'valgrind.xml')))
            removeTempDir()

if __name__ == '__main__':
    unittest.main()
