import mafIndelDistribution as MID
import numpy
import os
import re
import shutil
import sys
import unittest
myBinDir = os.path.normpath(os.path.dirname(sys.argv[0]))

class GenericObject:
    pass

class VerifyMafLineExtraction(unittest.TestCase):
    knownValues = [('s test1.chrA 0 10 + 100 ACTG---ACTGAC',
                    MID.MafLine('test1', 'chrA', 0, 10, 100, 1, 'ACTG---ACTGAC', 9)),
                   ('s test2.chrA+chrB.chr1.chr2+chr10 0 10 + 100 ACTG---ACTGAC',
                    MID.MafLine('test2', 'chrA+chrB.chr1.chr2+chr10', 0, 10, 100, 1, 'ACTG---ACTGAC', 9)),
                   ('s test2.chrA 0 17 + 100 -ACGTUMRWSYKVHDBXN-',
                    MID.MafLine('test2', 'chrA', 0, 17, 100, 1, '-ACGTUMRWSYKVHDBXN-', 9)),
                   ('s test1.chrA 0 10 - 100 ACTG---ACTGAC',
                    MID.MafLine('test1', 'chrA', 0, 10, 100, -1, 'ACTG---ACTGAC', 9)),
                   ]
    nameRegex = r'(.+?)\.(chr.+)'
    namePat = re.compile(nameRegex)
    options = GenericObject()
    options.genome1 = 'test1'
    options.genome2 = 'test2'
    options.speciesList = set(['test1', 'test2'])

    def test_extractMafLine(self):
        """ ensure that the function extractMafLine() parses lines as expected.
        """
        for pre, post in self.knownValues:
            x = MID.extractMafLine(self.namePat, pre, 9, self.options)
            if x is None or x == 'notOurGenome':
                self.assertEqual(x, post)
                continue
            #print ' '
            for v in ['genome', 'chrom', 'start', 'seqLength','totalLength',
                      'strand', 'sequence', 'lineno']:
                
                #print 'v:%11s pre:%s post[v]:%15s x[v]:%s' % (v, pre, post.__dict__[v], x.__dict__[v])
                self.assertEqual(post.__dict__[v], x.__dict__[v])

class VerifyMafRead(unittest.TestCase):
    header = '##maf version=1\n\n'
    footer = '\n'
    knownValues = (('case 1',
                    'a score=0\n'
                    's test1.chrA  0 10 + 100 ATGCATGCAT\n'
                    's test2.chrA  0 10 + 100 ATGCATGCAT\n',
                    {'test1' : {'chrA' : {'test2' : [(0, 10)]}},
                     'test2' : {'chrA' : {'test1' : [(0, 10)]}},}
                    ),
                   # case 2
                   ('case 2',
                    'a score=0\n'
                    's test1.chrA  0 10 - 100 ATGCATGCAT\n'
                    's test2.chrA  0 10 + 100 ATGCATGCAT\n',
                    {'test1' : {'chrA' : {'test2' : [(90, 100)]}},
                     'test2' : {'chrA' : {'test1' : [(0, 10)]}},}
                    ),
                   # case 3
                   ('case 3',
                    'a score=0\n'
                    's test1.chrA 90 10 - 100 ATGCATGCAT\n'
                    's test2.chrA  0 10 + 100 ATGCATGCAT\n',
                    {'test1' : {'chrA' : {'test2' : [(0, 10)]}},
                     'test2' : {'chrA' : {'test1' : [(0, 10)]}},}
                    ),
                   # case 4
                   ('case 4',
                    'a score=0\n'
                    's test1.chrA 90 10 - 100 ATGCATGCAT\n'
                    's test1.chrA 10 10 + 100 ATGCATGCAT\n'
                    's test2.chrA  0 10 + 100 ATGCATGCAT\n',
                    {'test1' : {'chrA' : {'test2' : [(0, 10), (10, 20)]}},
                     'test2' : {'chrA' : {'test1' : [(0, 10), (0, 10)]}},}
                    ),
                   # case 5
                   ('case 5',
                    'a score=0\n'
                    's test1.chrA 90 10 - 100 ATGCATGCAT\n'
                    's test1.chrA  9 10 + 100 ATGCATGCAT\n'
                    's test2.chrA 20 10 + 100 ATGCATGCAT\n'
                    's test2.chrA  0 10 + 100 ATGCATGCAT\n',
                    {'test1' : {'chrA' : {'test2' : [(9, 19), (0, 10), (9, 19), (0, 10)]}},
                     'test2' : {'chrA' : {'test1' : [(20, 30), (0, 10), (20, 30), (0, 10),]}},}
                    ),
                   # case 6
                   ('case 6',
                    'a score=0\n'
                    's test1.chrA 30 10 + 100 AT-GC-AT-GC-AT\n'
                    's test2.chrA  9 10 + 100 AT-GC-AT-GC-AT\n',
                    {'test1' : {'chrA' : {'test2' : [(30, 40)]}},
                     'test2' : {'chrA' : {'test1' : [(9, 19)]}},}
                    ),
                   # case 7
                   ('case 7',
                    'a score=0\n'
                    's test1.chrA 30  3 + 100 ATC----\n'
                    's test2.chrA  9  3 + 100 ----ACT\n',
                    {'test1' : {'chrA' : {'test2' : [(0, 0)]}},
                     'test2' : {'chrA' : {'test1' : [(0, 0)]}},}
                    ),
                   # case 8
                   ('case 8',
                    'a score=0\n'
                    's test1.chrA 30  3 + 100 ----ACT\n'
                    's test2.chrA  9  3 + 100 ACT----\n',
                    {'test1' : {'chrA' : {'test2' : [(0, 0)]}},
                     'test2' : {'chrA' : {'test1' : [(0, 0)]}},}
                    ),
                   # case 9
                   ('case 9',
                    'a score=0\n'
                    's test1.chrA  0  4 + 100 ACGT--\n'
                    's test2.chrA  0  4 + 100 --ACGT\n',
                    {'test1' : {'chrA' : {'test2' : [(2, 4)]}},
                     'test2' : {'chrA' : {'test1' : [(0, 2)]}},}
                    ),
                   # case 10
                   ('case 10',
                    'a score=0\n'
                    's test1.chrA  0  4 + 100 ACGT---\n'
                    's test2.chrA  0  4 + 100 ---ACGT\n',
                    {'test1' : {'chrA' : {'test2' : [(3, 4)]}},
                     'test2' : {'chrA' : {'test1' : [(0, 1)]}},}
                    ),
                   # case 11
                   ('case 11',
                    'a score=0\n'
                    's test1.chrA  0  4 - 100 ACGT--\n'
                    's test2.chrA  0  4 + 100 --ACGT\n',
                    {'test1' : {'chrA' : {'test2' : [(96, 98)]}},
                     'test2' : {'chrA' : {'test1' : [(0, 2)]}},}
                    ),
                   # case 12
                   ('case 12',
                    'a score=0\n'
                    's test2.chrA  0  4 + 100 --ACGT\n'
                    's test1.chrA  0  4 + 100 ACGT--\n',
                    {'test2' : {'chrA' : {'test1' : [(0, 2)]}},
                     'test1' : {'chrA' : {'test2' : [(2, 4)]}},}
                    ),
                   # case 13
                   ('case 13',
                    'a score=0\n'
                    's test2.chrA  0  4 + 100 ---ACGT\n'
                    's test1.chrA  0  4 + 100 ACGT---\n',
                    {'test2' : {'chrA' : {'test1' : [(0, 1)]}},
                     'test1' : {'chrA' : {'test2' : [(3, 4)]}},}
                    ),
                   # case 14
                   ('case 14',
                    'a score=0\n'
                    's test2.chrA  0  4 + 100 --ACGT\n'
                    's test1.chrA  0  4 - 100 ACGT--\n',
                    {'test2' : {'chrA' : {'test1' : [(0, 2)]}},
                     'test1' : {'chrA' : {'test2' : [(96, 98)]}},}
                    ),
                   # case 15
                   ('case 15',
                    'a score=0\n'
                    's test1.chrA  0  4 + 100 -ACGT--\n'
                    's test2.chrA  0  7 + 100 ACGTACG\n',
                    {'test1' : {'chrA' : {'test2' : [(0, 4)]}},
                     'test2' : {'chrA' : {'test1' : [(1, 5)]}},}
                    ),
                   # case 16
                   ('case 16',
                    'a score=0\n'
                    's test1.chrA  0  4 + 100 ACGT---\n'
                    's test2.chrA  0  7 + 100 ACGTACG\n',
                    {'test1' : {'chrA' : {'test2' : [(0, 4)]}},
                     'test2' : {'chrA' : {'test1' : [(0, 4)]}},}
                    ),
                   # case 17
                   ('case 17',
                    'a score=0\n'
                    's test1.chrA  0  4 + 100 ---ACGT\n'
                    's test2.chrA  0  7 + 100 ACGTACG\n',
                    {'test1' : {'chrA' : {'test2' : [(0, 4)]}},
                     'test2' : {'chrA' : {'test1' : [(3, 7)]}},}
                    ),
                   # case 18
                   ('case 18',
                    'a score=0\n'
                    's test1.chrA  0  7 + 100 ACGTACG\n'
                    's test2.chrA  0  4 + 100 -ACGT--\n',
                    {'test1' : {'chrA' : {'test2' : [(1, 5)]}},
                     'test2' : {'chrA' : {'test1' : [(0, 4)]}},}
                    ),
                   # case 19
                   ('case 19',
                    'a score=0\n'
                    's test1.chrA  0  7 + 100 ACGTACG\n'
                    's test2.chrA  0  4 + 100 ACGT---\n',
                    {'test1' : {'chrA' : {'test2' : [(0, 4)]}},
                     'test2' : {'chrA' : {'test1' : [(0, 4)]}},}
                    ),
                   # case 20
                   ('case 20',
                    'a score=0\n'
                    's test1.chrA  0  7 + 100 ACGTACG\n'
                    's test2.chrA  0  4 + 100 ---ACGT\n',
                    {'test1' : {'chrA' : {'test2' : [(3, 7)]}},
                     'test2' : {'chrA' : {'test1' : [(0, 4)]}},}
                    ),
                   # case 21
                   ('case 21',
                    'a score=0\n'
                    's test1.chrA  0 14 + 100 ACGTACGACGTACG\n'
                    's test2.chrA  0  8 + 100 ---ACGT---ACGT\n',
                    {'test1' : {'chrA' : {'test2' : [(3, 7), (10, 14)]}},
                     'test2' : {'chrA' : {'test1' : [(0, 8)]}},}
                    ),
                   # case 22
                   ('case 22',
                    'a score=0\n'
                    's test1.chrA  1 14 - 100 ACGTACGACGTACG\n'
                    's test2.chrA  1  8 - 100 ---ACGT---ACGT\n',
                    {'test1' : {'chrA' : {'test2' : [(85, 89), (92, 96)]}},
                     'test2' : {'chrA' : {'test1' : [(91, 99)]}},}
                    ),
                   # case 23
                   ('case 23',
                    'a score=0\n'
                    's test1.chrA  1 14 - 100 ACGTACGACGTACG\n'
                    's test2.chrA  1  8 - 100 ---ACGT---ACGT\n'
                    's test1.chrB  1 14 - 100 ACGTACGACGTACG\n',
                    {'test1' : {'chrA' : {'test2' : [(85, 89), (92, 96)]},
                                'chrB' : {'test2' : [(85, 89), (92, 96)]}},
                     'test2' : {'chrA' : {'test1' : [(91, 99), (91, 99)]}},
                     },
                    ),
                   )
    options = GenericObject()
    options.genome1 = 'test1'
    options.genome2 = 'test2'
    options.speciesList = set(['test1', 'test2'])
    if not os.path.exists('tempTestFiles'):
        os.mkdir('tempTestFiles')
    options.maf = os.path.join('tempTestFiles', 'testmaf.maf')
    
    def buildTruth(self, post):
        truth = {}
        for genomeA in post:
            truth[genomeA] = {}
            for chromA in post[genomeA]:
                truth[genomeA][chromA] = {}
                for genomeB in post[genomeA][chromA]:
                    truth[genomeA][chromA][genomeB] = numpy.zeros(100, dtype = numpy.uint16)
                    for start, stop in post[genomeA][chromA][genomeB]:
                        truth[genomeA][chromA][genomeB][start : stop] += 1
        return truth
    
    def test_oneWay(self):
        """ readMaf() should parse a maf file as expected
        """ 
        for name, pre, post in self.knownValues:
            f = open(self.options.maf, 'w')
            f.write('%s%s%s' % (self.header, pre, self.footer))
            f.close()
            
            alignments = MID.readMaf(self.options.maf, self.options)
            trueAlignments = self.buildTruth(post)
            
            # print name
            for g1 in alignments:
               self.assertTrue(g1 in trueAlignments)
               for c1 in alignments[g1]:
                   self.assertTrue(c1 in trueAlignments[g1])
                   for g2 in alignments[g1][c1]:
                       self.assertTrue(g2 in trueAlignments[g1][c1])
                       # print alignments[g1][c1][g2][c2]
                       # print trueAlignments[g1][c1][g2][c2]
                       # print sum(alignments[g1][c1][g2][c2] == trueAlignments[g1][c1][g2][c2])
                       self.assertTrue((sum(alignments[g1][c1][g2] == trueAlignments[g1][c1][g2]) == 
                                        len(alignments[g1][c1][g2])))
        shutil.rmtree(os.path.dirname(self.options.maf))

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
                predicted += MID.analyzeOne(a, self.options)
            self.assertEqual(predicted, gapTruth)
            
    def test_oneWayNoEdges(self):
        """ analyze() should return the correct values for a given vector when NoEdges is turned on
        """ 
        for name, arrayList, gapTruth, coveredTruth in self.knownValuesNoEdges:
            predicted = []
            for a in arrayList:
                predicted += MID.analyzeOne(a, self.optionsNoEdges)
            self.assertEqual(predicted, gapTruth)

    def test_coverage(self):
        """ calcBasesCovered() should return the correct values for a given vector
        """
        for name, arrayList, gapTruth, coveredTruth in self.knownValues:
            predicted = []
            for a in arrayList:
                predicted.append(MID.calcBasesCovered(a, self.options))
            self.assertEqual(predicted, coveredTruth)

if __name__ == '__main__':
   unittest.main()
