import unittest
import os

from sonLib.bioio import system
from sonLib.bioio import getTempFile

class TestCase(unittest.TestCase):
    def testRepeatBed(self):
        tempFile = getTempFile(rootDir=os.getcwd())
        tempFile2 = getTempFile(rootDir=os.getcwd())
        fileHandle = open(tempFile, 'w')
        fileHandle.write(">hello boo\nacTGACCCCgtcgAAcAAccc\n>foo\nAaaAAAAAAA")
        fileHandle.close()
        system("eval_getRepeatBed %s %s" % (tempFile, tempFile2))
        fileHandle = open(tempFile2, 'r')
        fn = lambda (i, j, k) : (i, int(j), int(k))
        j = [ fn(i.split()) for i in fileHandle.readlines() ]
        print j
        assert j == [ ("hello", 0, 2), ("hello", 9, 13), ("hello", 15, 16), ("hello", 18, 21), ("foo", 1, 3) ]
        os.remove(tempFile)
        os.remove(tempFile2)
        
if __name__ == '__main__':
    unittest.main()