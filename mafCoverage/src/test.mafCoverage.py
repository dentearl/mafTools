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

class CoverageTest(unittest.TestCase):
  def testCoverage_0(self):
    """ mafCoverage should output numbers that accurately reflect the
    coverage between sequences.
    """
    mtt.makeTempDirParent()
    tmpDir = os.path.abspath(mtt.makeTempDir('testCoverage_0'))
    parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    cmd = mtt.genericValgrind(tmpDir)
    cmd.append(os.path.abspath(os.path.join(parent, 'test', 'mafCoverage')))
    outpipes = [os.path.join('/dev', 'null')]
    mtt.recordCommands([cmd], tmpDir, outPipes=outpipes)
    mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
    self.assertTrue(False)
    mtt.removeDir(tmpDir)

  def _testMemory0(self):
    """ If valgrind is installed on the system, check for memory
    related errors (0).
    """
    valgrind = mtt.which('valgrind')
    if valgrind is None:
      return
    mtt.makeTempDirParent()
    tmpDir = os.path.abspath(mtt.makeTempDir('memory_0'))
    parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    cmd = mtt.genericValgrind(tmpDir)
    cmd.append(os.path.abspath(os.path.join(parent, 'test', 'mafCoverage')))
    outpipes = [os.path.join('/dev', 'null')]
    mtt.recordCommands([cmd], tmpDir, outPipes=outpipes)
    mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
    self.assertTrue(mtt.noMemoryErrors(os.path.join(tmpDir, 'valgrind.xml')))
    self.assertTrue(False)
    mtt.removeDir(tmpDir)


class CuTestMemory(unittest.TestCase):
  def test_CuTestMemory(self):
    """ If valgrind is installed on the system, check for memory
    related errors in CuTests.
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
