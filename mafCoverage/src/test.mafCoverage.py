##################################################
# Copyright (C) 2013 by
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
import math
import numpy as np
import os
import random
import sys
import unittest
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]),
                                             '../../include/')))
import mafToolsTest as mtt
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]),
                                             '../../mafValidator/src/')))
import mafValidator as mafval

g_headers = ['''##maf version=1 scoring=tba.v8
# tba.v8 (((human chimp) baboon) (mouse rat))

''',
             '''##maf version=1 scoring=tba.v8
# tba.v8 (((human chimp) baboon) (mouse rat))
''']

# known good : tuples (target species name, maf sequence,
#   dictionary key: species value: coverage)
g_knownGood = [('target', '''a score=0.0 status=test.input
s target.chr1  0 10 + 10 ACGTACGTAC
s seqB.chr1    0 10 + 10 ACGTACGTAC
s seqC.chr1    0  5 + 10 A----CGT-C

a score=0.0 status=test.input
s seqA.chr2  0 10 + 10 ACGTACGTAC
s seqB.chr3  0 10 - 10 ACGTACGTAC
s seqC.chr4  0  5 + 10 A----CGT-C

''',
                {'seqA': 0.0,
                 'seqB': 1.0,
                 'seqC': 0.5,
                 'target': 0.0,}),
               ('target',
                '''a score=0.0 status=test.input
s target.chr1  0 10 + 10 ACGTACGTAC
s seqB.chr1    0 10 + 10 ACGTACGTAC
s seqC.chr1    0  8 + 10 A-GT-CCGTC
s seqA.chr2    0  2 + 10 A--------C

a score=0.0 status=test.input
s seqA.chr2  0 10 + 10 ACGTACGTAC
s seqB.chr3  0 10 - 10 ACGTACGTAC
s seqC.chr4  0 10 + 10 ACGTACGTAC

''',
                {'seqA': 0.2,
                 'seqB': 1.0,
                 'seqC': 0.8,
                 'target': 0.0,}),
               ('target',
                '''a score=0.0 status=test.input
s target.chr1  0 10 + 10 ACGTACGTAC
s seqB.chr1    0 10 + 10 ACGTACGTAC
s seqC.chr1    0  8 + 10 A-GT-CCGTC
s seqA.chr2    0  2 + 10 A--------C

a score=0.0 status=test.input
s seqA.chr2  0 10 + 10 ACGTACGTAC
s seqB.chr3  0 10 - 10 ACGTACGTAC
s seqC.chr4  0 10 + 10 ACGTACGTAC

a score=0.0 status=test.input
s target.chr2  0 10 + 10 ACGTACGTAC
s seqD.chr5    0 10 + 10 ACGTACGTAC
s seqA.chr5    0  2 + 10 A--------C

''',
                {'seqA': 0.2,
                 'seqB': 0.5,
                 'seqC': 0.4,
                 'seqD': 0.5,
                 'target': 0.0,}),
               ('target',
                '''a score=0.0 status=test.input
s target.chr1  0 9 +  9 ACGTACGTA
s seqB.chr1    0 9 + 10 ACGTACGTA
s seqC.chr1    0 6 + 10 A--TAC-TA
s seqA.chr1    0 3 + 10 A--T----A

a score=0.0 status=test.input
s seqA.chr2  0 10 + 10 ACGTACGTAC
s seqB.chr3  0 10 - 10 ACGTACGTAC
s seqZ.chr4  0 10 + 10 ACGTACGTAC

''',
                {'seqA': 1.0 / 3.0,
                 'seqB': 1.0,
                 'seqC': 2.0 / 3.0,
                 'seqZ': 0.0,
                 'target': 0.0,}),
              ]

def createPredictionDict(path):
  f = open(path, 'r')
  predDict = {}
  # file format columns:
  # Target Species - querySpecies - length of Q genome - coverage - nCoverage
  for line in f:
    line = line.strip()
    if line.startswith('#') or line.startswith('targetSpecies'):
      continue
    tabs = line.split()
    predDict[tabs[1]] = float(tabs[3])
  return predDict


def AlmostEqual(a, b):
  ''' returns true of a and b are within tolerance of one another
  '''
  tolerance = 0.00001
  return abs(a - b) < tolerance


def dictionariesAreIdentical(dict_a, dict_b):
  ''' check that these dictionaries are identical
  '''
  for a in dict_a:
    if a not in dict_b:
      print 'a not in dict_b: %s, ' % (a, str(dict_b))
      return False
    if not AlmostEqual(dict_a[a], dict_b[a]):
      print 'dict_a[%s] %f != dict_b[%s] %f' % (a, dict_a[a], a, dict_b[a])
      return False
  for b in dict_b:
    if b not in dict_a:
      print 'b not in dict_a: %s, %s' % (b, str(dict_a))
      return False
    if not AlmostEqual(dict_b[b], dict_a[b]):
      print 'dict_b[%s] %f != dict_a[%s] %f' % (b, dict_b[b], b, dict_a[b])
      return False
  return True


def coverageIsCorrect(predictionPath, solutionDict):
  ''' read the predictionPath file, create a prediction dictionary,
  check that solution dictionary and prediction dictionary are identical
  '''
  predictionDict = createPredictionDict(predictionPath)
  if dictionariesAreIdentical(predictionDict, solutionDict):
    return True
  return False


class CoverageTest(unittest.TestCase):
  def testCoverage_0(self):
    """ mafCoverage should output numbers that accurately reflect the coverage between sequences.
    """
    mtt.makeTempDirParent()
    tmpDir = os.path.abspath(mtt.makeTempDir('testCoverage_0'))
    parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    customOpts = mafval.GenericValidationOptions()
    for target, mafSeq, solutionDict in g_knownGood:
      testMaf = mtt.testFile(os.path.join(os.path.abspath(tmpDir), 'test.maf'),
                             mafSeq, g_headers)
      cmd = [os.path.abspath(os.path.join(parent, 'test', 'mafCoverage'))]
      cmd += ['--maf', os.path.abspath(os.path.join(tmpDir, 'test.maf')),
              '--species', target]
      outpipes = [os.path.join(tmpDir, 'output.txt')]
      mtt.recordCommands([cmd], tmpDir, outPipes=outpipes)
      mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
      self.assertTrue(mafval.validateMaf(os.path.join(tmpDir, 'test.maf'),
                                         customOpts))
      self.assertTrue(coverageIsCorrect(os.path.join(tmpDir, 'output.txt'),
                                        solutionDict))
    mtt.removeDir(tmpDir)

  def testMemory0(self):
    """ If valgrind is installed on the system, check for memory related errors (0).
    """
    valgrind = mtt.which('valgrind')
    if valgrind is None:
      return
    mtt.makeTempDirParent()
    tmpDir = os.path.abspath(mtt.makeTempDir('memory_0'))
    parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    customOpts = mafval.GenericValidationOptions()
    for target, mafSeq, solutionDict in g_knownGood:
      testMaf = mtt.testFile(os.path.join(os.path.abspath(tmpDir), 'test.maf'),
                             mafSeq, g_headers)
      cmd = mtt.genericValgrind(tmpDir)
      cmd.append(os.path.abspath(os.path.join(parent, 'test', 'mafCoverage')))
      cmd += ['--maf', os.path.abspath(os.path.join(tmpDir, 'test.maf')),
              '--species', target]
      outpipes = [os.path.join('/dev', 'null')]
      mtt.recordCommands([cmd], tmpDir, outPipes=outpipes)
      mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
      self.assertTrue(mafval.validateMaf(os.path.join(tmpDir, 'test.maf'),
                                         customOpts))
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
