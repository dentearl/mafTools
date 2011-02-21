import unittest
import os

"""
Copyright (C) 2009-2011 by 
Dent Earl (dearl@soe.ucsc.edu, dentearl@gmail.com)
Benedict Paten (benedict@soe.ucsc.edu, benedictpaten@gmail.com)
Mark Diekhans (markd@soe.ucsc.edu)
... and other members of the Reconstruction Team of David Haussler's 
lab (BME Dept. UCSC).

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

"""

from sonLib.bioio import system
from sonLib.bioio import getTempFile

class TestCase(unittest.TestCase):
    def testRepeatBed(self):
        tempFile = getTempFile(rootDir=os.getcwd())
        tempFile2 = getTempFile(rootDir=os.getcwd())
        fileHandle = open(tempFile, 'w')
        fileHandle.write(">hello boo\nacTGACCCCgtcgAAcAAccc\n>foo\nAaaAAAAAAA")
        fileHandle.close()
        system("getRepeatBed %s %s" % (tempFile, tempFile2))
        fileHandle = open(tempFile2, 'r')
        fn = lambda (i, j, k) : (i, int(j), int(k))
        j = [ fn(i.split()) for i in fileHandle.readlines() ]
        print j
        assert j == [ ("hello", 0, 2), ("hello", 9, 13), ("hello", 15, 16), ("hello", 18, 21), ("foo", 1, 3) ]
        os.remove(tempFile)
        os.remove(tempFile2)
        
if __name__ == '__main__':
    unittest.main()
