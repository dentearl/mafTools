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
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), '../../lib/')))
import mafToolsTest as mtt

g_targetSeq1 = 'target.chr0'
g_targetSeq2 = 'target2.chr0'

g_headers = ['''##maf version=1 scoring=tba.v8
# tba.v8 (((human chimp) baboon) (mouse rat))

''',
             '''##maf version=1 scoring=tba.v8
# tba.v8 (((human chimp) baboon) (mouse rat))
''',
             '''##maf version=1

''',]
g_bedString = '''target.chr0 10 20
target.chr0 25 30
target.chr0 62 68
target.chr0 90 95
'''
g_overlappingBlocks = [('''a score=0
s target.chr0        38 13 +       100 gcagctgaaaaca
s target2.chr0       38 13 +       200 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA
s name.chr1           0 10 +       100 ATGT---ATGCCG
s name2.chr1         50 10 +       100 ATGT---ATGCCG
s name3.chr9         50 10 +       100 ATGTA---TGCCG
s name4.chr&         50 10 +       100 ATG---TATGCCG
s name5              50 10 +       100 ATGT---ATGCCG

''', 13, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
            ('''a score=0
s name                0 10 +       100 ATGTATGC---CG
s target2.chr0        0  1 -       200 A------------
s name2.chr1         50 10 +       100 ATGTATG---CCG
s name3.chr9         50 10 +       100 ATGTATG---CCG
s name4.chr&         50 10 +       100 ATGTATG---CCG
s name5              50 10 +       100 ATGTATG---CCG
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s target.chr0         5 10 -       100 ATGTATG---CCG

''', 1, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 1, 0, 0, 0, 0, 0]),
            ('''a score=0
s name               10 10 +       100 ATGTAT---GCCG
s name2.chr1         50 10 +       100 ATGTAT---GCCG
s name3.chr9         50 10 +       100 ATGTAT---GCCG
s name4.chr&         50 10 +       100 ATGTAT---GCCG
s name5              50 10 +       100 ATGTAT---GCCG
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s target.chr0        62  9 +       100 gca---gaa-aca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s hg16.chr7    27707221 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s target2.chr0       62  9 +       200 gcagaa----aca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA

''', 6, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 1, 1, 1, 0, 0, 0, 1, 1,
         1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0,]),
            ]
g_coverageLines = set(['target.chr01006220', 'target2.chr02002320'])
g_coverageLinesWild = set(['target.*2010010220', 'target2.chr02002320', 'target.chr01006220', 'target.chr320000400',
                           ])
g_coverageLinesWildBed = set(['#Overall', '#SequenceSourceLengthAlignedPos.Coverage',
                              'target.*2010010220', 'target2.chr02002320',
                              '#IndividualSequences', 'target.chr01006220',
                              'target.chr320000400', 'target2.chr02002320',
                              '#BEDRegions', '#SequenceRegionLengthInRegionsOutofRegionsCoverage',
                              'target.*26416', 'target.chr026416'])
g_nonOverlappingBlocks = [('''a score=23262.0
s hg18.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon         116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG

''', 0, None),
                          ('''a score=23262.0
s hg18.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon         116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
s target.chr3      1234 40 +     20000 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG

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
                          ('''a degree=4
s simChimp.chrA     15393341 348 + 53121445 ATATTGAGGAGCAGGATGGGTATAGAGGCCCTGACCCATTAATGTGTAAGCACTAGGCAGCTGGGAGATACCCCAGAGGGCGGGGTCACTGAATTCACTGGCCCACCACTGTAAATACATTCTAACCAGTGGGTTTAGGGCTCTGTGCATTAGAACCACTCTGAAGAAGTGTAACACACCACCTAGTGAGCTGCCGGGCCGCCAGCAACTTCTTTTTCCCACATGACCCATGCAAGCCCGTGATTTCTCCCTGGTACATGATATTTGGGATTCCAGGGACCTAATGGAGCATGCTATTCCTGTGTTAGTTATCACTTCGAAGGGGGTGCAAGAGTGTAAGTAATGGGT
s simGorilla.chrA   15595743 348 + 53120926 ATATTGAGGAGCAGGATGGGTATAGAAGCCCTGACCCATTAATGTGTAAGCACTAGGCAGCTGGGAGATACTCCAGAGGGAGGGGTCACTGAATTCACTGGCCCACCACTGTAAATACATTCTAACCAGTGGGTTTAGGGCTCGGTGCATTAGAACCACCCTGAAGAAGTGTAACGCACCACCTAGTGAGCTGCCGGGCCGCCAGCAACTTCTTTTTCCCACATGACCCATGCATGCCCGTGATTTCTCCCTGGTACATGGTTTTTGGGATTCCAGGGACCTAATGGAGCATACTATTCCTGTGTTAGTTATCACTTCGAAGGGGGTGCGAGAGTGTAAGTAATGGGT
s simHuman.chrA     36713600 348 - 53106993 ATATTGAGGAGCAGGATGGGTATAGAAGCCCTGACCTATTAATGTGTAAGCACTAGGCAGCTGGGCGATACCCCAGAGGGAGGGGTCACTGAATTCACTGGCCCACCACTGTAAATACATTCTAACCAGTGGGTTTAGGGCTCTGTGCATTAGAACCACCCTGAAGAAGAGTAACGCACCACCTAGTGAGCTGCCGGGCCGCCAGCAAGTTCTTTTTCCCACATGACCCATGCAAGCCCGTGATTTCTCCCTGGTACATGATATTTGGGATTCCAGGGACCTAATGGAGCATGCTATTCCTGTGTTAGTTATCACTTCGAAGGGGGTGCAAGAGTGTAAGTAATGGGT
s simOrang.chrE       126390 348 + 37692687 ATTTTGAGGAGCAGGATGGGTATAGAAGCCCTGACCCATTAATGTGTGAGCTCTAGGCAGCTTGGAGATACTGCAGAGGGAGGGGTCACTGAATTCACTGGCCCACCACTGTAAATACATACTAACCGGTGGGTTTAGGGCTCTGTGCATTAGAACCACCCTGAGGAAGTGTAACGCACCACCTAGTGAGCTGCCGGGCCACCAGCAACTTCTTTTTCCCACATGACCCATGCAAGCCCGTGATTTCTCCCTGGTACATGATCTTTGGGATTCCAGGGACCTAATGGCGGATGCTATTCCTGTGTTAGTTATCACTTCGAAGGGGGCGCAAGAGTGTAAGTAATGGGT
s target.chr0     0 30 - 100 AAAAAAAAAAAAAAAAAAAAAAAAAAAAA------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T

''', 0, [0] * 100),
    ]

def coverageIsCorrect(filename, lineSet):
  f = open(filename, 'r')
  for line in f:
    line = line.strip()
    if line.startswith('#'):
      if line == '# Bins':
        # This function does not test bin values
        return True
      continue
    if not line:
      continue
    line = line.split()
    if len(line):
      line.pop() # throw away the first coverage number
      if len(line) == 5:
        line.pop() # throw away the second coverage number
    line = ''.join(line)
    if line not in lineSet:
      print '%s not in lineSet' % line
      return False
  return True


def binningIsCorrect(filename, true_array):
  """given an outputfile name and the expected array, verify the output
  is correct
  """
  f = open(filename, 'r')
  is_bin_section = False
  obs_array = []
  for line in f:
    line = line.strip()
    if line.startswith('#'):
      if line == '# Bins':
        is_bin_section = True
      continue
    if not line:
      continue
    if not is_bin_section:
      continue
    data = line.split()
    assert(len(data) == 2)
    obs_array.append(float(data[1]))
  obs_array = np.array(obs_array)
  if len(obs_array) != len(true_array):
    print 'len(obs_array) %d != len(true_array) %d' % (len(obs_array),
                                                       len(true_array))
    print 'obs:  ', obs_array
    print 'true: ', true_array
    return False
  if not np.allclose(obs_array, true_array):
    for i in xrange(0, len(obs_array)):
      if obs_array[i] != true_array[i]:
        print 'obs_array[%d] %e != true_array[%d] %e' % (i, obs_array[i],
                                                         i, true_array[i])
        print 'obs:  ', obs_array
        print 'true: ', true_array
        return False
  return True


def compressList(a_list, bin_length):
  """Takes a binary (0s and 1s) list and then reduces the length
  to be ceil(len(a_list) / bin_length) and sums all the values from
  a_list to go in the new list's larger buckets.
  """
  b_list = np.zeros(math.ceil(len(a_list) / float(bin_length)))
  cur_bin = -1
  for i in xrange(0, len(a_list)):
    if not (i % bin_length):
      cur_bin += 1
    b_list[cur_bin] +=  a_list[i]
  b_list /= bin_length
  return b_list

class CoverageTest(unittest.TestCase):
  def testCoverage_0(self):
    """ mafPairCoverage should output numbers that accurately reflect the coverage between sequences.
    """
    mtt.makeTempDirParent()
    for i in xrange(0, 10):
      shuffledBlocks = []
      tmpDir = os.path.abspath(mtt.makeTempDir('coverage_0'))
      order = [1] * len(g_overlappingBlocks) + [0] * len(g_nonOverlappingBlocks)
      random.shuffle(order)
      random.shuffle(g_overlappingBlocks)
      random.shuffle(g_nonOverlappingBlocks)
      j, k = 0, 0
      total = 0
      for isOverlapping in order:
        if isOverlapping:
          shuffledBlocks.append(g_overlappingBlocks[j][0])
          total += g_overlappingBlocks[j][1]
          j += 1
        else:
          shuffledBlocks.append(g_nonOverlappingBlocks[k][0])
          k += 1
      testMaf = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test_%d.maf' % i)),
                             ''.join(shuffledBlocks), g_headers)
      parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
      cmd = [os.path.abspath(os.path.join(parent, 'test', 'mafPairCoverage'))]
      cmd += ['--maf', os.path.abspath(os.path.join(tmpDir, 'test_%d.maf' % i)),
              '--seq1', g_targetSeq1, '--seq2', g_targetSeq2,
              ]
      outpipes = [os.path.abspath(os.path.join(tmpDir, 'coverage.txt'))]
      mtt.recordCommands([cmd], tmpDir, outPipes=outpipes)
      mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
      self.assertTrue(
        coverageIsCorrect(os.path.join(tmpDir, 'coverage.txt'), g_coverageLines))
      mtt.removeDir(tmpDir)

  def testCoverage_1(self):
    """ mafPairCoverage should output numbers that accurately reflect the coverage between sequences.
    """
    mtt.makeTempDirParent()
    for i in xrange(0, 10):
      shuffledBlocks = []
      tmpDir = os.path.abspath(mtt.makeTempDir('coverage_1'))
      order = [1] * len(g_overlappingBlocks) + [0] * len(g_nonOverlappingBlocks)
      random.shuffle(order)
      random.shuffle(g_overlappingBlocks)
      random.shuffle(g_nonOverlappingBlocks)
      j, k = 0, 0
      total = 0
      for isOverlapping in order:
        if isOverlapping:
          shuffledBlocks.append(g_overlappingBlocks[j][0])
          total += g_overlappingBlocks[j][1]
          j += 1
        else:
          shuffledBlocks.append(g_nonOverlappingBlocks[k][0])
          k += 1
      testMaf = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test_%d.maf' % i)),
                             ''.join(shuffledBlocks), g_headers)
      parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
      cmd = [os.path.abspath(os.path.join(parent, 'test', 'mafPairCoverage'))]
      cmd += ['--maf', os.path.abspath(os.path.join(tmpDir, 'test_%d.maf' % i)),
              '--seq1', 'target.*', '--seq2', g_targetSeq2,
              ]
      outpipes = [os.path.abspath(os.path.join(tmpDir, 'coverage.txt'))]
      mtt.recordCommands([cmd], tmpDir, outPipes=outpipes)
      mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
      self.assertTrue(coverageIsCorrect(os.path.join(tmpDir, 'coverage.txt'), g_coverageLinesWild))
      mtt.removeDir(tmpDir)

  def testCoverageBed_0(self):
    """ mafPairCoverage should be able to get the correct output given bed-based region intervals
    """
    mtt.makeTempDirParent()
    for i in xrange(0, 10):
      shuffledBlocks = []
      tmpDir = os.path.abspath(mtt.makeTempDir('bed_0'))
      order = [1] * len(g_overlappingBlocks) + [0] * len(g_nonOverlappingBlocks)
      random.shuffle(order)
      random.shuffle(g_overlappingBlocks)
      random.shuffle(g_nonOverlappingBlocks)
      j, k = 0, 0
      total = 0
      for isOverlapping in order:
        if isOverlapping:
          shuffledBlocks.append(g_overlappingBlocks[j][0])
          total += g_overlappingBlocks[j][1]
          j += 1
        else:
          shuffledBlocks.append(g_nonOverlappingBlocks[k][0])
          k += 1
      testMaf = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test_%d.maf' % i)),
                             ''.join(shuffledBlocks), g_headers)
      testBed = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'intervals.bed')), g_bedString)
      parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
      cmd = [os.path.abspath(os.path.join(parent, 'test', 'mafPairCoverage'))]
      cmd += ['--maf', os.path.abspath(os.path.join(tmpDir, 'test_%d.maf' % i)),
              '--seq1', 'target.*', '--seq2', g_targetSeq2, '--bed',
              os.path.abspath(os.path.join(tmpDir, 'intervals.bed'))
              ]
      outpipes = [os.path.abspath(os.path.join(tmpDir, 'coverage.txt'))]
      mtt.recordCommands([cmd], tmpDir, outPipes=outpipes)
      mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
      self.assertTrue(coverageIsCorrect(os.path.join(tmpDir, 'coverage.txt'), g_coverageLinesWildBed))
      mtt.removeDir(tmpDir)

  def testCoverageBinning_0(self):
    """ mafPairCoverage should output bin values that accurately reflect the coverage between sequences.
    """
    mtt.makeTempDirParent()
    for i in xrange(0, 10):
      shuffledBlocks = []
      true_array = np.zeros(100)
      tmpDir = os.path.abspath(mtt.makeTempDir('binning_0'))
      order = [1] * len(g_overlappingBlocks) + [0] * len(g_nonOverlappingBlocks)
      random.shuffle(order)
      random.shuffle(g_overlappingBlocks)
      random.shuffle(g_nonOverlappingBlocks)
      j, k = 0, 0
      total = 0
      for isOverlapping in order:
        if isOverlapping:
          shuffledBlocks.append(g_overlappingBlocks[j][0])
          total += g_overlappingBlocks[j][1]
          true_array = np.add(true_array, g_overlappingBlocks[j][2])
          j += 1
        else:
          shuffledBlocks.append(g_nonOverlappingBlocks[k][0])
          k += 1
      testMaf = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test_%d.maf' % i)),
                             ''.join(shuffledBlocks), g_headers)
      parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
      cmd = [os.path.abspath(os.path.join(parent, 'test', 'mafPairCoverage'))]
      cmd += ['--maf', os.path.abspath(os.path.join(tmpDir, 'test_%d.maf' % i)),
              '--seq1', g_targetSeq1, '--seq2', g_targetSeq2, '--bin_start=0',
              '--bin_end=99', '--bin_length=1',
              ]
      outpipes = [os.path.abspath(os.path.join(tmpDir, 'coverage.txt'))]
      mtt.recordCommands([cmd], tmpDir, outPipes=outpipes)
      mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
      self.assertTrue(
        coverageIsCorrect(os.path.join(tmpDir, 'coverage.txt'), g_coverageLines))
      self.assertTrue(
        binningIsCorrect(os.path.join(tmpDir, 'coverage.txt'), true_array))
      mtt.removeDir(tmpDir)

  def testCoverageBinning_1(self):
    """ mafPairCoverage should output bin values that accurately reflect the coverage between sequences.
    """
    mtt.makeTempDirParent()
    for i in xrange(0, 10):
      shuffledBlocks = []
      true_array = np.zeros(100)
      tmpDir = os.path.abspath(mtt.makeTempDir('binning_0'))
      order = [1] * len(g_overlappingBlocks) + [0] * len(g_nonOverlappingBlocks)
      random.shuffle(order)
      random.shuffle(g_overlappingBlocks)
      random.shuffle(g_nonOverlappingBlocks)
      j, k = 0, 0
      total = 0
      for isOverlapping in order:
        if isOverlapping:
          shuffledBlocks.append(g_overlappingBlocks[j][0])
          total += g_overlappingBlocks[j][1]
          true_array = np.add(true_array, g_overlappingBlocks[j][2])
          j += 1
        else:
          shuffledBlocks.append(g_nonOverlappingBlocks[k][0])
          k += 1
      testMaf = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test_%d.maf' % i)),
                             ''.join(shuffledBlocks), g_headers)
      parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
      cmd = [os.path.abspath(os.path.join(parent, 'test', 'mafPairCoverage'))]
      cmd += ['--maf', os.path.abspath(os.path.join(tmpDir, 'test_%d.maf' % i)),
              '--seq1', g_targetSeq1, '--seq2', g_targetSeq2, '--bin_start=0',
              '--bin_end=99', '--bin_length=2',
              ]
      outpipes = [os.path.abspath(os.path.join(tmpDir, 'coverage.txt'))]
      mtt.recordCommands([cmd], tmpDir, outPipes=outpipes)
      mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
      self.assertTrue(
        coverageIsCorrect(os.path.join(tmpDir, 'coverage.txt'), g_coverageLines))
      self.assertTrue(
        binningIsCorrect(os.path.join(tmpDir, 'coverage.txt'),
                         compressList(true_array, 2)))
      mtt.removeDir(tmpDir)

  def testCoverageBinning_2(self):
    """ mafPairCoverage should output bin values that accurately reflect the coverage between sequences, TEST A HUGE NUMBER OF BINS.
    """
    mtt.makeTempDirParent()
    for j in xrange(1, 100):  # bin lengths
      for i in xrange(0, 10):
        shuffledBlocks = []
        true_array = np.zeros(100)
        tmpDir = os.path.abspath(mtt.makeTempDir('binning_0'))
        order = [1] * len(g_overlappingBlocks) + [0] * len(g_nonOverlappingBlocks)
        random.shuffle(order)
        random.shuffle(g_overlappingBlocks)
        random.shuffle(g_nonOverlappingBlocks)
        j, k = 0, 0
        total = 0
        for isOverlapping in order:
          if isOverlapping:
            shuffledBlocks.append(g_overlappingBlocks[j][0])
            total += g_overlappingBlocks[j][1]
            true_array = np.add(true_array, g_overlappingBlocks[j][2])
            j += 1
          else:
            shuffledBlocks.append(g_nonOverlappingBlocks[k][0])
            k += 1
        testMaf = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test_%d.maf' % i)),
                               ''.join(shuffledBlocks), g_headers)
        parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        cmd = [os.path.abspath(os.path.join(parent, 'test', 'mafPairCoverage'))]
        cmd += ['--maf', os.path.abspath(os.path.join(tmpDir, 'test_%d.maf' % i)),
                '--seq1', g_targetSeq1, '--seq2', g_targetSeq2, '--bin_start=0',
                '--bin_end=99', '--bin_length=%d' % j,
                ]
        outpipes = [os.path.abspath(os.path.join(tmpDir, 'coverage.txt'))]
        mtt.recordCommands([cmd], tmpDir, outPipes=outpipes)
        mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
        self.assertTrue(
          coverageIsCorrect(os.path.join(tmpDir, 'coverage.txt'), g_coverageLines))
        self.assertTrue(
          binningIsCorrect(os.path.join(tmpDir, 'coverage.txt'),
                           compressList(true_array, j)))
        mtt.removeDir(tmpDir)

  def _testMemory0(self):
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
          shuffledBlocks.append(g_overlappingBlocks[j][0])
          j += 1
        else:
          shuffledBlocks.append(g_nonOverlappingBlocks[k][0])
          k += 1
      testMaf = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test_%d.maf' % i)),
                             ''.join(shuffledBlocks), g_headers)
      parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
      cmd = mtt.genericValgrind(tmpDir)
      cmd.append(os.path.abspath(os.path.join(parent, 'test', 'mafPairCoverage')))
      cmd += ['--maf', os.path.abspath(os.path.join(tmpDir, 'test_%d.maf' % i)),
              '--seq1', g_targetSeq1, '--seq2', g_targetSeq2,
              ]
      outpipes = [os.path.join('/dev', 'null')]
      mtt.recordCommands([cmd], tmpDir, outPipes=outpipes)
      mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
      self.assertTrue(mtt.noMemoryErrors(os.path.join(tmpDir, 'valgrind.xml')))
      mtt.removeDir(tmpDir)

  def _testMemory1(self):
    """ If valgrind is installed on the system, check for memory related errors (0).
    """
    valgrind = mtt.which('valgrind')
    if valgrind is None:
      return
    mtt.makeTempDirParent()
    for i in xrange(0, 10):
      shuffledBlocks = []
      tmpDir = os.path.abspath(mtt.makeTempDir('memory_1'))
      order = [1] * len(g_overlappingBlocks) + [0] * len(g_nonOverlappingBlocks)
      random.shuffle(order)
      random.shuffle(g_overlappingBlocks)
      random.shuffle(g_nonOverlappingBlocks)
      j, k = 0, 0
      for isOverlapping in order:
        if isOverlapping:
          shuffledBlocks.append(g_overlappingBlocks[j][0])
          j += 1
        else:
          shuffledBlocks.append(g_nonOverlappingBlocks[k][0])
          k += 1
      testMaf = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'test_%d.maf' % i)),
                             ''.join(shuffledBlocks), g_headers)
      parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
      cmd = mtt.genericValgrind(tmpDir)
      cmd.append(os.path.abspath(os.path.join(parent, 'test', 'mafPairCoverage')))
      cmd += ['--maf', os.path.abspath(os.path.join(tmpDir, 'test_%d.maf' % i)),
              '--seq1', 'target.*', '--seq2', g_targetSeq2,
              ]
      outpipes = [os.path.join('/dev', 'null')]
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
