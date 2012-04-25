import mafCoveragePickleAnalysis as MCPA
import numpy
import os
import re
import shutil
import sys
import unittest
myBinDir = os.path.normpath(os.path.dirname(sys.argv[0]))

class GenericObject:
    pass

class VerifyAnalysis(unittest.TestCase):
    # gaps defined to exist anywhere, even at the edges of alignments
    knownValues = [('case 1',
                    [numpy.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1], dtype = numpy.uint16),
                     numpy.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1], dtype = numpy.uint16)], 
                    [], [10, 10] ),
                   ('case 2',
                   [numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype = numpy.uint16),
                     numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype = numpy.uint16)],
                    [10, 10], [0, 0], ),
                   ('case 3',
                   [numpy.array([0, 0, 0, 0, 1, 1, 1, 0, 0, 0], dtype = numpy.uint16),
                    numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype = numpy.uint16)],
                    [4, 3, 10], [3, 0], ),
                   ('case 4',
                   [numpy.array([0, 0, 0, 0, 1, 1, 1, 0, 0, 0], dtype = numpy.uint16),
                    numpy.array([0, 1, 0, 1, 0, 1, 0, 1, 0, 0], dtype = numpy.uint16)],
                    [4, 3, 1, 1, 1, 1, 2], [3, 4], ),
                   ('case 5',
                   [numpy.array([0, 0, 0, 0, 1, 1, 1, 0, 0, 0], dtype = numpy.uint16),
                    numpy.array([1, 0, 0, 0, 1, 1, 1, 1, 0, 0], dtype = numpy.uint16),
                    numpy.array([0, 1, 0, 1, 0, 1, 0, 1, 0, 0], dtype = numpy.uint16)],
                    [4, 3, 3, 2, 1, 1, 1, 1, 2], [3, 5, 4], ),
                   ]
    options = GenericObject()
    options.noEdges = False

    # gaps defined to exist only between regions of alignment
    knownValuesNoEdges = [('case 1',
                           [numpy.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1], dtype = numpy.uint16),
                            numpy.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1], dtype = numpy.uint16)], 
                           [], [10, 10] ),
                          ('case 2',
                           [numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype = numpy.uint16),
                            numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype = numpy.uint16)],
                           [], [0, 0], ),
                          ('case 3',
                           [numpy.array([0, 0, 0, 0, 1, 1, 1, 0, 0, 0], dtype = numpy.uint16),
                            numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype = numpy.uint16)],
                           [], [3, 0], ),
                          ('case 4',
                           [numpy.array([0, 0, 0, 0, 1, 1, 1, 0, 0, 0], dtype = numpy.uint16),
                            numpy.array([0, 1, 0, 1, 0, 1, 0, 1, 0, 0], dtype = numpy.uint16)],
                           [1, 1, 1], [3, 4], ),
                          ('case 5',
                           [numpy.array([0, 0, 0, 0, 1, 1, 1, 0, 0, 0], dtype = numpy.uint16),
                            numpy.array([1, 0, 0, 0, 1, 1, 1, 1, 0, 0], dtype = numpy.uint16),
                            numpy.array([0, 1, 0, 1, 0, 1, 0, 1, 0, 0], dtype = numpy.uint16)],
                           [3, 1, 1, 1], [3, 5, 4], ),
                          ]
    optionsNoEdges = GenericObject()
    optionsNoEdges.noEdges = True
    
    def buildArray(intervalList):
        a = numpy.zeros(100, dtype = numpy.uint16)
        for start, stop in intervalList:
            a[start : stop] += 1
        return a
    
    def test_oneWay(self):
        """ analyze() should return the correct values for a given vector
        """ 
        for name, arrayList, gapTruth, coveredTruth in self.knownValues:
            predicted = []
            for a in arrayList:
                predicted += MCPA.analyzeOne(a, self.options)
            self.assertEqual(predicted, gapTruth)
            
    def test_oneWayNoEdges(self):
        """ analyze() should return the correct values for a given vector when NoEdges is turned on
        """ 
        for name, arrayList, gapTruth, coveredTruth in self.knownValuesNoEdges:
            predicted = []
            for a in arrayList:
                predicted += MCPA.analyzeOne(a, self.optionsNoEdges)
            self.assertEqual(predicted, gapTruth)

    def test_coverage(self):
        """ calcBasesCovered() should return the correct values for a given vector
        """
        for name, arrayList, gapTruth, coveredTruth in self.knownValues:
            predicted = []
            for a in arrayList:
                predicted.append(MCPA.calcBasesCovered(a, self.options))
            self.assertEqual(predicted, coveredTruth)

if __name__ == '__main__':
   unittest.main()
