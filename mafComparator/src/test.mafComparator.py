##############################
# Copyright (C) 2009-2012 by
# Dent Earl (dearl@soe.ucsc.edu, dentearl@gmail.com)
# Benedict Paten (benedict@soe.ucsc.edu, benedictpaten@gmail.com)
# Mark Diekhans (markd@soe.ucsc.edu)
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
##############################
import xml.etree.ElementTree as ET
import os
import shutil
import subprocess
import unittest
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), '../../include/')))
import mafToolsTest as mtt

g_headers = ['''##maf version=1
''',
             '''##maf version=1 scoring=tba.v8
# tba.v8 (((human chimp) baboon) (mouse rat))
''',
             '''##maf version=1 scoring=tba.v8
# tba.v8 (((human chimp) baboon) (mouse rat))

''',]

def xmlBedRegionPassed(filename, totalTrue, totalTrueInInterval):
    tree = ET.parse(filename)
    homTests = tree.findall('homologyTests')
    if totalTrue != int(homTests[0].find('aggregateResults').find('all').attrib['totalTrue']):
        return False
    if totalTrueInInterval is None:
        if homTests[0].find('aggregateResults').find('A') is not None:
            return False
    else:
        if totalTrueInInterval != int(homTests[0].find('aggregateResults').find('A').attrib['totalTrue']):
            return False
    return True
def getAggregateResult(filename, name):
    tree = ET.parse(filename)
    homTests = tree.findall('homologyTests')
    return int(homTests[0].find('aggregateResults').find('all').attrib[name])

class KnownValuesTest(unittest.TestCase):
    # knownValues contains quad-tuples,
    # maf1, maf2, totalTrue (comparing maf1 as fileA to maf2), 
    # totalFalse (comparing maf1 as fileA to maf2)
    knownValues = [('''a score=0
s A 0 10 + 20 ACGTACGTAC
s B 0 10 + 20 ATGTACGTAC

''', '''a score=0
s A 0 10 + 20 ACGTACGTAC
s B 0 10 + 20 ATGTACGTAC

''', 10, 0),
                   ('''a score=0
s A 0 10 + 20 ACGTACGTAC
s B 0 10 + 20 ATGTACGTAC

''', '''a score=0
s A 10 10 + 20 GTACGTACGT
s B 10 10 + 20 ATGTACGTAC

''', 0, 10),
                   ('''a score=0
s A 0 10 + 10 ACGTACGTAC
s B 0 10 + 10 ATGTACGTAC

''', '''a score=0
s A 0 10 - 10 GTACGTACGT
s B 0 10 - 10 ATGTACGTAC

''', 10, 0),
                   ('''a score=0
s A 0 5 + 10 ACGTA
s B 0 5 + 10 ATGTA

''', '''a score=0
s A  5 5 - 10 TACGT
s B  5 5 - 10 TACAT

''', 5, 0),
                   ('''a score=0
s A 0 10 + 20 ACGTACGTAC
s B 0 10 + 20 ATGTACGTAC

''', '''a score=0
s A 10 10 - 20 GTACGTACGT
s B 10 10 - 20 GTACGTACAT

''', 10, 0),
                   ('''a score=0
s A 8 10 + 20 ACGTACGTAC
s B 8 10 + 20 ATGTACGTAC

''', '''a score=0
s A 0 10 + 20 ACGTACGTAC
s B 0 10 + 20 ATGTACGTAC

''', 2, 8),
                   ]
    def test_knownValues(self):
        """ mafComparator should return correct results for hand-calculable problems
        """
        mtt.makeTempDirParent()
        tmpDir = os.path.abspath(mtt.makeTempDir('knownValues'))
        for maf1, maf2, totalTrue, totalFalse in self.knownValues:
            testMaf1 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf1.maf')), 
                                    maf1, g_headers)
            testMaf2 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf2.maf')), 
                                    maf2, g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = [os.path.abspath(os.path.join(parent, 'test', 'mafComparator')),
                   '--maf1', os.path.abspath(os.path.join(tmpDir, 'maf1.maf')),
                   '--maf2', os.path.abspath(os.path.join(tmpDir, 'maf2.maf')),
                   '--out', os.path.abspath(os.path.join(tmpDir, 'output.xml')),
                   '--samples=1000', '--logLevel=critical',
                   ]
            mtt.recordCommands([cmd], tmpDir)
            mtt.runCommandsS([cmd], tmpDir)
            passedTT = totalTrue == getAggregateResult(os.path.abspath(os.path.join(tmpDir, 'output.xml')), 'totalTrue')
            passedTF = totalFalse == getAggregateResult(os.path.abspath(os.path.join(tmpDir, 'output.xml')), 'totalFalse')
            self.assertTrue(passedTT and passedTF)
        mtt.removeDir(tmpDir)
    def test_memory_2(self):
        """ mafComparator should be memory clean for known values
        """
        valgrind = mtt.which('valgrind')
        if valgrind is None:
            return
        mtt.makeTempDirParent()        
        tmpDir = os.path.abspath(mtt.makeTempDir('knownValues'))
        for maf1, maf2, totalTrue, totalFalse in self.knownValues:
            testMaf1 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf1.maf')), 
                                    maf1, g_headers)
            testMaf2 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf2.maf')), 
                                    maf2, g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = mtt.genericValgrind(tmpDir)
            cmd += [os.path.abspath(os.path.join(parent, 'test', 'mafComparator')),
                    '--maf1', os.path.abspath(os.path.join(tmpDir, 'maf1.maf')),
                    '--maf2', os.path.abspath(os.path.join(tmpDir, 'maf2.maf')),
                    '--out', os.path.abspath(os.path.join(tmpDir, 'output.xml')),
                    '--samples=1000', '--logLevel=critical',
                    ]
            mtt.recordCommands([cmd], tmpDir)
            mtt.runCommandsS([cmd], tmpDir)
            passed = mtt.noMemoryErrors(os.path.join(tmpDir, 'valgrind.xml'))
            self.assertTrue(passed)
        mtt.removeDir(tmpDir)

class BedParsingTest(unittest.TestCase):
    # knownValues contains quad-tuples,
    # maf1, maf2, bed, threshold
                   # test 1 - 0 % bed coverage
    knownValues = [('a score=0\n'
                    's test1.chr0 0 20 + 100 ACGTACGTACGTACGTACGT\n'
                    's test2.chr0 0 20 + 100 ACGTACGTACGTACGTACGT\n\n',
                    'a score=0\n'
                    's test1.chr0 0 20 + 100 ACGTACGTACGTACGTACGT\n'
                    's test2.chr0 0 20 + 100 ACGTACGTACGTACGTACGT\n\n',
                    'test1.chr0\t30\t100\n',
                    20, 0),
                   # test 2 - no bed file
                   ('a score=0\n'
                    's test1.chr0 0 20 + 100 ACGTACGTACGTACGTACGT\n'
                    's test2.chr0 0 20 + 100 ACGTACGTACGTACGTACGT\n\n',
                    'a score=0\n'
                    's test1.chr0 0 20 + 100 ACGTACGTACGTACGTACGT\n'
                    's test2.chr0 0 20 + 100 ACGTACGTACGTACGTACGT\n\n',
                    '',
                    20, None),
                   # test 3 - 50% bed coverage
                   ('a score=0\n'
                    's test1.chr0 0 20 + 100 ACGTACGTACGTACGTACGT\n'
                    's test2.chr0 0 20 + 100 ACGTACGTACGTACGTACGT\n\n',
                    'a score=0\n'
                    's test1.chr0 0 20 + 100 ACGTACGTACGTACGTACGT\n'
                    's test2.chr0 0 20 + 100 ACGTACGTACGTACGTACGT\n\n',
                    'test1.chr0\t11\t100\n',
                    20, 9),
                   # test 3 - 100% bed coverage
                   ('a score=0\n'
                    's test1.chr0 0 20 + 100 ACGTACGTACGTACGTACGT\n'
                    's test2.chr0 0 20 + 100 ACGTACGTACGTACGTACGT\n\n',
                    'a score=0\n'
                    's test1.chr0 0 20 + 100 ACGTACGTACGTACGTACGT\n'
                    's test2.chr0 0 20 + 100 ACGTACGTACGTACGTACGT\n\n',
                    'test1.chr0\t1\t100\n',
                    20, 19),
                   ]
    def test_bedParsing(self):
        """ mafComparator should parse a bed file and use the intervals for testing
        """
        mtt.makeTempDirParent()
        tmpDir = os.path.abspath(mtt.makeTempDir('bedParsing'))
        for maf1, maf2, bed, totalTrue, totalTrueInInterval in self.knownValues:
            testMaf1 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf1.maf')), 
                                    maf1, g_headers)
            testMaf2 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf2.maf')), 
                                    maf2, g_headers)
            testBed = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'bed.bed')), 
                                   bed, [''])
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = [os.path.abspath(os.path.join(parent, 'test', 'mafComparator')),
                   '--maf1', os.path.abspath(os.path.join(tmpDir, 'maf1.maf')),
                   '--maf2', os.path.abspath(os.path.join(tmpDir, 'maf2.maf')),
                   '--out', os.path.abspath(os.path.join(tmpDir, 'output.xml')),
                   '--samples=1000', '--logLevel=critical',
                   ]
            if bed != '':
                cmd += ['--bedFiles', os.path.abspath(os.path.join(tmpDir, 'bed.bed'))]
            mtt.recordCommands([cmd], tmpDir)
            mtt.runCommandsS([cmd], tmpDir)
            passed = xmlBedRegionPassed(os.path.abspath(os.path.join(tmpDir, 'output.xml')), 
                                        totalTrue, totalTrueInInterval)
            self.assertTrue(passed)
        mtt.removeDir(tmpDir)
    def test_memory_0(self):
        """ mafComparator should be memory clean for bed parsing examples
        """
        valgrind = mtt.which('valgrind')
        if valgrind is None:
            return
        mtt.makeTempDirParent()
        tmpDir = os.path.abspath(mtt.makeTempDir('memory_0'))
        for maf1, maf2, bed, totalTrue, totalTrueInInterval in self.knownValues:
            testMaf1 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf1.maf')), 
                                    maf1, g_headers)
            testMaf2 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf2.maf')), 
                                    maf2, g_headers)
            testBed = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'bed.bed')), 
                                   bed, [''])
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = mtt.genericValgrind(tmpDir)
            cmd += [os.path.abspath(os.path.join(parent, 'test', 'mafComparator')),
                    '--maf1', os.path.abspath(os.path.join(tmpDir, 'maf1.maf')),
                    '--maf2', os.path.abspath(os.path.join(tmpDir, 'maf2.maf')),
                    '--out', os.path.abspath(os.path.join(tmpDir, 'output.xml')),
                    '--samples=1000', '--logLevel=critical',
                    ]
            if bed != '':
                cmd += ['--bedFiles', os.path.abspath(os.path.join(tmpDir, 'bed.bed'))]
            mtt.recordCommands([cmd], tmpDir)
            mtt.runCommandsS([cmd], tmpDir)
            passed = mtt.noMemoryErrors(os.path.join(tmpDir, 'valgrind.xml'))
            self.assertTrue(passed)
        mtt.removeDir(tmpDir)
class randomSeedTests(unittest.TestCase):
    knownValues = [('a score=0\n'
                    's test1.chr0 0 20 + 100 ACGTACGTACGTACGTACGT\n'
                    's test2.chr0 0 20 + 100 ACGTACGTACGTACGTACGT\n',
                    'a score=0\n'
                    's test1.chr0 0 20 + 100 ACGTACGTACGTACGTACGT\n'
                    's test2.chr0 0 20 + 100 ACGTACGTACGTACGTACGT\n',
                    ),
                   ('a score=0\n'
                    's test1.chr0 0 20 + 100 ACGTACGTACGTACGTACGT\n'
                    's test2.chr0 0 20 + 100 ACGTACGTACGTACGTACGT\n',
                    'a score=0\n'
                    's test1.chr0 10 30 + 100 ACGTACGTACGTACGTACGT\n'
                    's test2.chr0 10 30 + 100 ACGTACGTACGTACGTACGT\n',
                    ),
                   ('a score=0\n'
                    's test1.chr0 0 20 + 100 ACGTACGTACGTACGTACGT\n'
                    's test2.chr0 0 20 + 100 ACGTACGTACGTACGTACGT\n',
                    'a score=0\n'
                    's test1.chr0 40 60 + 100 ACGTACGTACGTACGTACGT\n'
                    's test2.chr0 40 60 + 100 ACGTACGTACGTACGTACGT\n',
                    ),
                   ]
    def test_seedTesting(self):
        """ mafComparator should have replicable runs via the --seed command
        """
        mtt.makeTempDirParent()
        tmpDir = os.path.abspath(mtt.makeTempDir('seedTesting'))
        for maf1, maf2  in self.knownValues:
            testMaf1 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf1.maf')), 
                                    maf1, g_headers)
            testMaf2 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf2.maf')), 
                                    maf2, g_headers)
            cmd = ['mafComparator', 
                   '--mafFile1', os.path.abspath(os.path.join(tmpDir, 'maf1.maf')), 
                   '--mafFile2', os.path.abspath(os.path.join(tmpDir, 'maf2.maf')),
                   '--outputFile', os.path.join(tmpDir, 'output.xml'),
                   '--sampleNumber=10', '--logLevel=critical']
            mtt.recordCommands([cmd], tmpDir)
            mtt.runCommandsS([cmd], tmpDir)
            tree = ET.parse(os.path.join(tmpDir, 'output.xml'))
            ac = tree.getroot()
            seed = int(ac.attrib['seed'])
            cmd.append('--seed=%d' % seed)
            origHomTests = tree.findall('homologyTests')
            for i in xrange(0, 10):
                mtt.recordCommands([cmd], tmpDir)
                mtt.runCommandsS([cmd], tmpDir)
                tree = ET.parse(os.path.join(tmpDir, 'output.xml'))
                ac = tree.getroot()
                homTests = tree.findall('homologyTests')
                self.assertEqual(seed, int(ac.attrib['seed']))
                for elm in ['totalTrue', 'totalFalse', 'average']:
                    self.assertEqual(homTests[0].find('aggregateResults').find('all').attrib[elm],
                                     origHomTests[0].find('aggregateResults').find('all').attrib[elm])
                    self.assertEqual(homTests[1].find('aggregateResults').find('all').attrib[elm],
                                     origHomTests[1].find('aggregateResults').find('all').attrib[elm])
        mtt.removeDir(tmpDir)
    def test_memory_1(self):
        """ mafComparator should be memory clean for seed testing examples
        """
        valgrind = mtt.which('valgrind')
        if valgrind is None:
            return
        mtt.makeTempDirParent()
        tmpDir = os.path.abspath(mtt.makeTempDir('memory_1'))
        for maf1, maf2  in self.knownValues:
            testMaf1 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf1.maf')), 
                                    maf1, g_headers)
            testMaf2 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf2.maf')), 
                                    maf2, g_headers)
            cmd = mtt.genericValgrind(tmpDir)
            cmd += ['mafComparator', 
                    '--mafFile1', os.path.abspath(os.path.join(tmpDir, 'maf1.maf')), 
                    '--mafFile2', os.path.abspath(os.path.join(tmpDir, 'maf2.maf')),
                    '--outputFile', os.path.join(tmpDir, 'output.xml'),
                    '--sampleNumber=10', '--logLevel=critical']
            mtt.recordCommands([cmd], tmpDir)
            mtt.runCommandsS([cmd], tmpDir)
            self.assertTrue(mtt.noMemoryErrors(os.path.join(tmpDir, 'valgrind.xml')))
            for i in xrange(0, 4):
                mtt.recordCommands([cmd], tmpDir)
                mtt.runCommandsS([cmd], tmpDir)
                self.assertTrue(mtt.noMemoryErrors(os.path.join(tmpDir, 'valgrind.xml')))
        mtt.removeDir(tmpDir)
class CuTestTests(unittest.TestCase):
    def test_CuTestTests(self):
        """ Yo dawg, I heard you liked unit tests so I put some unit tests in your unit test so now you can unit test when you unit test.
        """
        mtt.makeTempDirParent()
        tmpDir = os.path.abspath(mtt.makeTempDir('CuTestAllTests'))
        cmd = mtt.genericValgrind(tmpDir)
        cmd += [os.path.abspath(os.path.join(os.curdir, 'test/allTests'))]
        outpipes = [os.path.join('/dev', 'null')]
        mtt.recordCommands([cmd], tmpDir)
        mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
        self.assertTrue(mtt.noMemoryErrors(os.path.join(tmpDir, 'valgrind.xml')))
        mtt.removeDir(tmpDir)

if __name__ == '__main__':
    unittest.main()
