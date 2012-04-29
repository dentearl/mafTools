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

def npa(a=None):
    """ numpy array, npa, shortcut function for creating a numpy array 
    with the correct datatype
    """
    if a is not None:
        if a.__class__ is int:
            a = [a]
    else:
        a = []
    return numpy.array(a, dtype=numpy.uint32)
def arraysEqual(a, b):
    if a.__class__ != b.__class__:
        # print 'classes differ, %s != %s' % (a.__class__, b.__class__)
        return False
    if a.__class__ == numpy.ndarray:
        if a.dtype != b.dtype:
            # print 'dtypes differ, %s != %s' % (a.dtype, b.dtype)
            return False
    for i in xrange(len(a)):
        if a[i] != b[i]:
            # print 'values differ, i %d: %d != %d' % (i, a[i], b[i])
            return False
    return True

class VerifyAnalysis(unittest.TestCase):
    """
    """
    # knownValues is a list of quads which are (name, list of input arrays, 
    # expected dict output, list of expected number of covered bases for each array
    knownValues = [('case 0',
                    [numpy.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1], dtype=numpy.uint16)],
                    {0: npa(), 1: npa(10), 'zerosAndOnes': npa(10), 'greaterThanOne': npa()}, 
                    [10]),
                   ('case 1',
                    [numpy.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1], dtype=numpy.uint16),
                     numpy.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1], dtype=numpy.uint16)], 
                    {0: npa(), 1: npa([10, 10]), 'zerosAndOnes': npa([10, 10]), 'greaterThanOne': npa()},
                   [10, 10],),
                   ('case 2',
                   [numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype=numpy.uint16),
                     numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype=numpy.uint16)],
                    {0: npa([10, 10]), 'zerosAndOnes': npa([10, 10]), 'greaterThanOne': npa()},
                    [0, 0],),
                   ('case 3',
                   [numpy.array([0, 0, 0, 0, 1, 1, 1, 0, 0, 0], dtype=numpy.uint16),
                    numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype=numpy.uint16)],
                    {0: npa([3, 4, 10]), 1: npa(3), 'zerosAndOnes': npa([10, 10]), 'greaterThanOne': npa()},
                    [3, 0], ),
                   ('case 4',
                   [numpy.array([0, 0, 0, 0, 1, 1, 1, 0, 0, 0], dtype=numpy.uint16),
                    numpy.array([0, 1, 0, 1, 0, 1, 0, 1, 0, 0], dtype=numpy.uint16)],
                    {0: npa([1, 1, 1, 1, 2, 3, 4]), 1: npa([1, 1, 1, 1, 3]), 
                     'zerosAndOnes': npa([10, 10]), 'greaterThanOne': npa()},
                    [3, 4], ),
                   ('case 5',
                   [numpy.array([0, 0, 0, 0, 1, 1, 1, 0, 0, 0], dtype=numpy.uint16),
                    numpy.array([1, 0, 0, 0, 1, 1, 1, 1, 0, 0], dtype=numpy.uint16),
                    numpy.array([0, 1, 0, 1, 0, 1, 0, 1, 0, 0], dtype=numpy.uint16)],
                    {0: npa([1, 1, 1, 1, 2, 2, 3, 3, 4]), 1: npa([1, 1, 1, 1, 1, 3, 4]), 
                     'zerosAndOnes': npa([10, 10, 10]), 'greaterThanOne': npa()},
                    [3, 5, 4], ),
                   ('case 6',
                   [numpy.array([0, 2, 2, 0, 1, 1, 1, 0, 3, 3], dtype=numpy.uint16),
                    numpy.array([1, 0, 0, 0, 1, 1, 1, 1, 0, 7], dtype=numpy.uint16),
                    numpy.array([0, 1, 0, 2, 2, 2, 2, 1, 0, 0], dtype=numpy.uint16)],
                    {0: npa([1, 1, 1, 1, 1, 1, 2, 3]), 1: npa([1, 1, 1, 3, 4]), 2: npa([2, 4]),
                     3: npa([2]), 4: npa(), 5: npa(), 6: npa(), 7: npa(1), 
                     'zerosAndOnes': npa([1, 3, 3, 5, 9]), 'greaterThanOne': npa([1, 2, 2, 4,])},
                    [7, 6, 6], ),
                   ]
    options = GenericObject()
    def test_oneWay(self):
        """ analyzeOne() should return the correct dictionary for a given number of alignment arrays
        """ 
        for name, arrayList, cnvTruth, coveredTruth in self.knownValues:
            cnvs = {}
            for a in arrayList:
                MCPA.analyzeOne(a, cnvs, self.options)
            for c in cnvs:
                self.assertTrue(c in cnvTruth)
                if cnvs[c].shape == (0,): 
                    ## the way to check if cnvs[c] is empty
                    self.assertTrue(cnvTruth[c].shape == (0,))
                else:
                    cnvs[c].sort()
                    self.assertTrue(arraysEqual(cnvs[c], cnvTruth[c]))
    def test_coverage(self):
        """ calcBasesCovered() should return the correct values for a given vector
        """
        for name, arrayList, cnvTruth, coveredTruth in self.knownValues:
            predicted = []
            for a in arrayList:
                predicted.append(MCPA.calcBasesCovered(a, self.options))
            self.assertEqual(predicted, coveredTruth)
if __name__ == '__main__':
   unittest.main()
