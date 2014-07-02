##############################
# Copyright (C) 2009-2013 by
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
import xml.parsers.expat
import os
import shutil
import subprocess
import unittest
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), '../../lib/')))
import mafToolsTest as mtt

g_headers = ['''##maf version=1
''',
             '''##maf version=1 scoring=tba.v8
# tba.v8 (((human chimp) baboon) (mouse rat))
''',
             '''##maf version=1 scoring=tba.v8
# tba.v8 (((human chimp) baboon) (mouse rat))

''',]

knownValues = [('''a score=0
# test 0.0
s A 0 10 + 20 ACGTACGTAC
s B 0 10 + 20 ATGTACGTAC

''', '''a score=0
# test 0.1
s A 0 10 + 20 ACGTACGTAC
s B 0 10 + 20 ATGTACGTAC

''', 10, 0), # test 0
               ('''a score=0
# test 1.0
s A 0 10 + 20 ACGTACGTAC
s B 0 10 + 20 ATGTACGTAC

''', '''a score=0
# test 1.1
s A 10 10 + 20 GTACGTACGT
s B 10 10 + 20 ATGTACGTAC

''', 0, 10), # test 1
               ('''a score=0
# test 2.0
s A 0 10 + 10 ACGTACGTAC
s B 0 10 + 10 ATGTACGTAC

''', '''a score=0
# test 2.1
s A 0 10 - 10 GTACGTACGT
s B 0 10 - 10 ATGTACGTAC

''', 10, 0), # test 2
               ('''a score=0
# test 3.0
s A 0 5 + 10 ACGTA
s B 0 5 + 10 ATGTA

''', '''a score=0
# test 3.1
s A  5 5 - 10 TACGT
s B  5 5 - 10 TACAT

''', 5, 0), # test 3
               ('''a score=0
# test 4.0
s A 0 10 + 20 ACGTACGTAC
s B 0 10 + 20 ATGTACGTAC

''', '''a score=0
# test 4.1
s A 10 10 - 20 GTACGTACGT
s B 10 10 - 20 GTACGTACAT

''', 10, 0), # test 4
               ('''a score=0
# test 5.0
s A 8 10 + 20 ACGTACGTAC
s B 8 10 + 20 ATGTACGTAC

''', '''a score=0
# test 5.1
s A 0 10 + 20 ACGTACGTAC
s B 0 10 + 20 ATGTACGTAC

''', 2, 8), # test 5
               ('''a score=0
# test 6.0
s A 0 10 + 20 ACGTACGTAC
s B 0 10 + 20 ATGTACGTAC
s C 0 10 + 20 ACGTACGTAC
s D 0 10 + 20 ATGTACGTAC

''', '''a score=0
# test 6.1
s A 0 10 + 20 ACGTACGTAC
s B 0 10 + 20 ATGTACGTAC
s C 0 10 + 20 ACGTACGTAC
s D 0 10 + 20 ATGTACGTAC

''', 60, 0), # test 6
               ('''a score=0
# test 7.0
s A 0 10 + 20 ACGTACGTAC
s B 0 10 + 20 ATGTACGTAC
s C 0 10 + 20 ACGTACGTAC
s D 0 10 + 20 ATGTACGTAC

''', '''a score=0
# test 7.1
s A  0 10 + 20 ACGTACGTAC
s B  0 10 + 20 ATGTACGTAC
s C 10 10 + 20 ACGTACGTAC
s D  0 10 + 20 ATGTACGTAC

''', 30, 30), # test 7
               ('''a score=0
# test 8.0
s A 0 10 + 20 ACGTACGTAC
s B 0  5 + 20 A-----GTAC
s C 0 10 + 20 ACGTACGTAC
s D 0 10 + 20 ATGTACGTAC

''', '''a score=0
# test 8.1
s A 0 10 + 20 ACGTACGTAC
s B 0 10 + 20 ATGTACGTAC
s C 0 10 + 20 ACGTACGTAC
s D 0 10 + 20 ATGTACGTAC

''',  33, 12), # test 8
               ('''a score=0
# test 9.0
s A 0 10 + 20 ACGTACGTAC
s B 0 10 + 20 ATGTACGTAC
s C 0 10 + 20 ACGTACGTAC
s D 0 10 + 20 ATGTACGTAC

''', '''a score=0
# test 9.1
s A 0 10 + 20 ACGTACGTAC
s B 0  5 + 20 A-----GTAC
s C 0 10 + 20 ACGTACGTAC
s D 0 10 + 20 ATGTACGTAC

''', 33, 27), # test 9
               ('''a score=0
# test 10.0
s A 0 10 + 20 ACGTACGTAC
s B 0 10 + 20 ATGTACGTAC
s C 0 10 + 20 ACGTACGTAC
s D 0 10 + 20 ATGTACGTAC

''', '''a score=0
# test 10.1
s A 0 10 + 20 ACGTACGTAC
s B 0 10 + 20 ATGTACGTAC
s C 0 10 + 20 ACGTACGTAC
s D 1 10 + 20 ATGTACGTAC

''', 30, 30), # test 10
               ('''a score=0
# test 11.0
s z 0 10 + 20 ACGTACGTAC
s A 0 10 + 20 ACGTACGTAC
s y 0 10 + 20 ACGTACGTAC
s x 0 10 + 20 ACGTACGTAC
s w 0 10 + 20 ACGTACGTAC
s B 0 10 + 20 ATGTACGTAC
s t 0 10 + 20 ACGTACGTAC
s C 0 10 + 20 ACGTACGTAC
s D 0 10 + 20 ATGTACGTAC

''', '''a score=0
# test 11.1
s k 0 10 + 20 ACGTACGTAC
s l 0 10 + 20 ACGTACGTAC
s A 0 10 + 20 ACGTACGTAC
s B 0 10 + 20 ATGTACGTAC
s m 0 10 + 20 ACGTACGTAC
s n 0 10 + 20 ACGTACGTAC
s C 0 10 + 20 ACGTACGTAC
s D 0 10 + 20 ATGTACGTAC
s p 0 10 + 20 ACGTACGTAC
s q 0 10 + 20 ACGTACGTAC

''', 60, 0), # test 11
               ('''a score=0
# test 12.0
s Anc3.0                       75484 13 + 607219 ACCAAGGAGTATG
s simMouse_chr6_simMouse.chr6 548643 13 - 636262 ACCACGGGGTATG

''', '''a score=0
# test 12.1
s Anc3.0                      531722 13 - 607219 CATACTCCTTGGT
s simMouse_chr6_simMouse.chr6  87606 13 + 636262 CATACCCCGTGGT

''', 13, 0), # test 12
               ('''a score=0
# test 13.1
s Anc3.0                      531722 13 - 607219 CATACTCCTTGGT
s simMouse_chr6_simMouse.chr6  87606 13 + 636262 CATACCCCGTGGT

''', '''a score=0
# test 13.0
s Anc3.0                       75484 13 + 607219 ACCAAGGAGTATG
s simMouse_chr6_simMouse.chr6 548643 13 - 636262 ACCACGGGGTATG

''', 13, 0), # test 13
               ]
knownValuesLegit = ['A:20,B:20', 'A:20,B:20', 'A:10,B:10', 'A:10,B:10', 'A:20,B:20',
                    'A:20,B:20', 'A:20,B:20,C:20,D:20', 'A:20,B:20,C:20,D:20', 'A:20,B:20,C:20,D:20', 
                    'A:20,B:20,C:20,D:20', 'A:20,B:20,C:20,D:20', 'A:20,B:20,C:20,D:20', 
                    'Anc3.0:607219,simMouse_chr6_simMouse.chr6:636262', 
                    'Anc3.0:607219,simMouse_chr6_simMouse.chr6:636262'
                    ]
knownValuesLegitFail = ['A:20,B:21', 'A:21,B:20', 'A:10,B:11', 'A:11,B:10', 'A:20,B:21',
                        'A:21,B:20', 'A:20,B:20,C:20,D:21', 'A:20,B:20,C:20,D:21', 'A:20,B:20,C:21,D:20', 
                        'A:20,B:21,C:20,D:20', 'A:21,B:20,C:20,D:20', 'A:20,B:20,C:20,D:21', 
                        'Anc3.0:607219,simMouse_chr6_simMouse.chr6:636263', 
                        'Anc3.0:607218,simMouse_chr6_simMouse.chr6:636262'
                        ]
knownValuesNumberOfPairs = ['10,10', '10,10', '10,10', '5,5', '10,10',
                            '10,10', '60,60', '60,60', '45,60', '60,45',
                            '60,60', '60,60', '13,13', '13,13'
                            ]
knownValuesNumberOfPairsFail = ['100,10', '100,10', '10,100', '5,50', '100,10',
                                '100,10', '60,600', '600,60', '450,60', '60,450',
                                '60,600', '600,60', '130,13', '13,130'
                                ]
knownValuesBed = [('a score=0\n'
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
knownValuesSeed = [('a score=0\n'
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
knownValuesNear = [('''a score=0
s A 0 10 + 20 ACGTACGTAC
s B 0 10 + 20 ATGTACGTAC
s C 0 10 + 20 ACGTACGTAC
s D 0 10 + 20 ATGTACGTAC

''', '''a score=0
s A 0 10 + 20 ACGTACGTAC
s B 0 10 + 20 ATGTACGTAC
s C 0 10 + 20 ACGTACGTAC
s D 0 10 + 20 ATGTACGTAC

''', 0, 60, 0),
                   ('''a score=0
s A 0 10 + 20 ACGTACGTAC
s B 0 10 + 20 ATGTACGTAC
s C 0 10 + 20 ACGTACGTAC
s D 0 10 + 20 ATGTACGTAC

''', '''a score=0
s A 0 10 + 20 ACGTACGTAC
s B 0 10 + 20 ATGTACGTAC
s C 0 10 + 20 ACGTACGTAC
s D 1 10 + 20 ATGTACGTAC

''', 1, 60, 0),
                   ('''a score=0
s A 0 10 + 20 ACGTACGTAC
s B 0 10 + 20 ATGTACGTAC
s C 0 10 + 20 ACGTACGTAC
s D 0 10 + 20 ATGTACGTAC

''', '''a score=0
s A 0 10 + 20 ACGTACGTAC
s B 0 10 + 20 ATGTACGTAC
s C 0 10 + 20 ACGTACGTAC
s D 1 10 + 20 ATGTACGTAC

''', 0, 30, 30),
                       ('''a score=0
s A 0 10 + 20 ACGTACGTAC
s B 0 10 + 20 ATGTACGTAC
s C 0 10 + 20 ACGTACGTAC
s D 0 10 + 20 ATGTACGTAC

''', '''a score=0
s A 0 10 + 20 ACGTACGTAC
s B 0 10 + 20 ATGTACGTAC
s C 0 10 + 20 ACGTACGTAC
s D 2 10 + 20 ATGTACGTAC

''', 2, 60, 0),
                   ('''a score=0
s A 0 10 + 20 ACGTACGTAC
s B 0 10 + 20 ATGTACGTAC
s C 0 10 + 20 ACGTACGTAC
s D 0 10 + 20 ATGTACGTAC

''', '''a score=0
s A 0 10 + 20 ACGTACGTAC
s B 0 10 + 20 ATGTACGTAC
s C 0 10 + 20 ACGTACGTAC
s D 2 10 + 20 ATGTACGTAC

''', 0, 30, 30),
                       ('''a score=0
s A 0 10 + 20 ACGTACGTAC
s B 0 10 + 20 ATGTACGTAC
s C 0 10 + 20 ACGTACGTAC
s D 0 10 + 20 ATGTACGTAC

''', '''a score=0
s A 0 10 + 20 ACGTACGTAC
s B 0 10 + 20 ATGTACGTAC
s C 0 10 + 20 ACGTACGTAC
s D 0 10 + 20 ATGTACGTAC

''', 5, 60, 0),
                       ]
knownValuesWiggles = [('''a score=0
s A 0 10 + 20 ACGTACGTAC
s B 0 10 + 20 ACGTACGTAC
''', 
                       '''a score=0
s A 0 10 + 20 ACGTACGTAC
s B 0 10 + 20 ACGTACGTAC
''', 
                       'A', 'B', 2,
                       [{'reference':'A',
                         'partner':'B',
                         'referenceLength':'20',
                         'numberOfBins':'10',
                         'binLength':'2',
                         'presentAtoB':'2,2,2,2,2,0,0,0,0,0',
                         'presentBtoA':'2,2,2,2,2,0,0,0,0,0',
                         'absentAtoB':'0,0,0,0,0,0,0,0,0,0',
                         'absentBtoA':'0,0,0,0,0,0,0,0,0,0',
                         },
                        ]),
                      # test 1
                      ('''a score=0
s A 0 10 + 20 ACGTACGTAC
s B 0 10 + 20 ACGTACGTAC
s C 0 10 + 20 ACGTACGTAC
''', 
                       '''a score=0
s A 0 10 + 20 ACGTACGTAC
s B 0 10 + 20 ACGTACGTAC
s C 0 10 + 20 ACGTACGTAC
''', 
                       'A', 'B', 2,
                       [{'reference':'A',
                         'partner':'B',
                         'referenceLength':'20',
                         'numberOfBins':'10',
                         'binLength':'2',
                         'presentAtoB':'2,2,2,2,2,0,0,0,0,0',
                         'presentBtoA':'2,2,2,2,2,0,0,0,0,0',
                         'absentAtoB':'0,0,0,0,0,0,0,0,0,0',
                         'absentBtoA':'0,0,0,0,0,0,0,0,0,0',
                         },
                        ]),
                      # test 2
                      ('''a score=0
s A.chr0 0 10 + 20 ACGTACGTAC
s B      0 10 + 20 ACGTACGTAC
''', 
                       '''a score=0
s A.chr0 0 10 + 20 ACGTACGTAC
s B      0 10 + 20 ACGTACGTAC
''', 
                       'A*', 'B', 2,
                       [{'reference':'A.chr0',
                         'partner':'B',
                         'referenceLength':'20',
                         'numberOfBins':'10',
                         'binLength':'2',
                         'presentAtoB':'2,2,2,2,2,0,0,0,0,0',
                         'presentBtoA':'2,2,2,2,2,0,0,0,0,0',
                         'absentAtoB':'0,0,0,0,0,0,0,0,0,0',
                         'absentBtoA':'0,0,0,0,0,0,0,0,0,0',
                         },
                        ]),
                      # test 3
                      ('''a score=0
s A.chr0 0 10 + 20 ACGTACGTAC
s B      0 10 + 20 ACGTACGTAC

s A.chr1 0 10 + 20 ACGTACGTAC
s B      0 10 + 20 ACGTACGTAC
''', 
                       '''a score=0
s A.chr0 0 10 + 20 ACGTACGTAC
s B      0 10 + 20 ACGTACGTAC

s A.chr1 0 10 + 20 ACGTACGTAC
s B      0 10 + 20 ACGTACGTAC
''', 
                       'A*', 'B', 2,
                       [{'reference':'A.chr0',
                         'partner':'B',
                         'referenceLength':'20',
                         'numberOfBins':'10',
                         'binLength':'2',
                         'presentAtoB':'2,2,2,2,2,0,0,0,0,0',
                         'presentBtoA':'2,2,2,2,2,0,0,0,0,0',
                         'absentAtoB':'0,0,0,0,0,0,0,0,0,0',
                         'absentBtoA':'0,0,0,0,0,0,0,0,0,0',
                         },
                        {'reference':'A.chr1',
                         'partner':'B',
                         'referenceLength':'20',
                         'numberOfBins':'10',
                         'binLength':'2',
                         'presentAtoB':'2,2,2,2,2,0,0,0,0,0',
                         'presentBtoA':'2,2,2,2,2,0,0,0,0,0',
                         'absentAtoB':'0,0,0,0,0,0,0,0,0,0',
                         'absentBtoA':'0,0,0,0,0,0,0,0,0,0',
                         },
                        ]),
                      # test 4
                      ('''a score=0
s A 0 10 + 20 ACGTACGTAC
s B 0 10 + 20 ACGTACGTAC
''', 
                       '''a score=0
s A 3 10 + 20 ACGTACGTAC
s B 3 10 + 20 ACGTACGTAC
''', 
                       'A', 'B', 2,
                       [{'reference':'A',
                         'partner':'B',
                         'referenceLength':'20',
                         'numberOfBins':'10',
                         'binLength':'2',
                         'presentAtoB':'0,1,2,2,2,0,0,0,0,0',
                         'presentBtoA':'0,1,2,2,2,0,0,0,0,0',
                         'absentAtoB' :'2,1,0,0,0,0,0,0,0,0',
                         'absentBtoA' :'0,0,0,0,0,2,1,0,0,0',
                         },
                        ]),
                      ]

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
def validXml(filename):
    try:
        tree = ET.parse(filename)
    except xml.parsers.expat.ExpatError:
        return False
    except ET.ParseError:
        return False
    return True
def validWiggleOutput(xml, valuesDictList):
    tree = ET.parse(xml)
    wiggles = tree.find('wigglePairs').findall('wigglePair')
    valid = True
    for wig in wiggles:
        matchFound = False
        for valuesDict in valuesDictList:
            if (valuesDict['reference'] == wig.attrib['reference'] and
                valuesDict['partner'] == wig.attrib['partner'] and
                valuesDict['referenceLength'] == wig.attrib['referenceLength'] and
                valuesDict['numberOfBins'] == wig.attrib['numberOfBins'] and
                valuesDict['binLength'] == wig.attrib['binLength'] and
                valuesDict['presentAtoB'] == wig.find('presentMaf1ToMaf2').text and
                valuesDict['presentBtoA'] == wig.find('presentMaf2ToMaf1').text and
                valuesDict['absentAtoB'] == wig.find('absentMaf1ToMaf2').text and
                valuesDict['absentBtoA'] == wig.find('absentMaf2ToMaf1').text):
                matchFound = True
        valid = valid and matchFound
    return valid
class XmlOutputValidation(unittest.TestCase):
    def test_validXmlOutput_(self):
        """ xml output file should be valid xml.
        """
        mtt.makeTempDirParent()
        tmpDir = os.path.abspath(mtt.makeTempDir('validateXmlOutput'))
        # test 1
        maf1, maf2, totalTrue, totalFalse = knownValues[0]
        testMaf1 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf1.maf')), 
                                maf1, g_headers)
        testMaf2 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf2.maf')), 
                                maf2, g_headers)
        parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        cmd = [os.path.abspath(os.path.join(parent, 'test', 'mafComparator')),
               '--maf1', os.path.abspath(os.path.join(tmpDir, 'maf1.maf')),
               '--maf2', os.path.abspath(os.path.join(tmpDir, 'maf2.maf')),
               '--out', os.path.abspath(os.path.join(tmpDir, 'outputVanilla.xml')),
               '--samples=1000', '--logLevel=critical',
               ]
        mtt.recordCommands([cmd], tmpDir)
        mtt.runCommandsS([cmd], tmpDir)
        passed = validXml(os.path.abspath(os.path.join(tmpDir, 'outputVanilla.xml')))
        self.assertTrue(passed)
        self.assertFalse(validXml(os.path.abspath(os.path.join(tmpDir, 'maf1.maf'))))
        # test 2
        maf1, maf2, bed, totalTrue, totalTrueInInterval = knownValuesBed[2]
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
               '--out', os.path.abspath(os.path.join(tmpDir, 'outputBed.xml')),
               '--samples=1000', '--logLevel=critical',
               ]
        cmd += ['--bedFiles', os.path.abspath(os.path.join(tmpDir, 'bed.bed'))]
        mtt.recordCommands([cmd], tmpDir)
        mtt.runCommandsS([cmd], tmpDir)
        passed = validXml(os.path.abspath(os.path.join(tmpDir, 'outputBed.xml')))
        # test 3
        maf1, maf2, totalTrue, totalFalse = knownValues[0]
        testMaf1 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf1.maf')), 
                                maf1, g_headers)
        testMaf2 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf2.maf')), 
                                maf2, g_headers)
        parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        cmd = [os.path.abspath(os.path.join(parent, 'test', 'mafComparator')),
               '--maf1', os.path.abspath(os.path.join(tmpDir, 'maf1.maf')),
               '--maf2', os.path.abspath(os.path.join(tmpDir, 'maf2.maf')),
               '--out', os.path.abspath(os.path.join(tmpDir, 'outputWiggles.xml')),
               '--samples=1000', '--logLevel=critical', '--wigglePairs=A:B', 'wiggleBinLength=2'
               ]
        mtt.recordCommands([cmd], tmpDir)
        mtt.runCommandsS([cmd], tmpDir)
        passed = validXml(os.path.abspath(os.path.join(tmpDir, 'outputWiggles.xml')))
        self.assertTrue(passed)
        mtt.removeDir(tmpDir)

class WigglePairsTest(unittest.TestCase):
    def test_wigglePairs(self):
        """ --wigglePairs options should produce expected and known output for given inputs
        """
        mtt.makeTempDirParent()
        tmpDir = os.path.abspath(mtt.makeTempDir('wigglePairs'))
        for maf1, maf2, ref, partner, binlen, valuesDictList in knownValuesWiggles:
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
                   '--wigglePairs=%s:%s' % (ref, partner),
                   '--wiggleBinLength=%d' % binlen,
                   ]
            mtt.recordCommands([cmd], tmpDir)
            mtt.runCommandsS([cmd], tmpDir)
            passed = validWiggleOutput(os.path.abspath(os.path.join(tmpDir, 'output.xml')), valuesDictList)
            self.assertTrue(passed)
        mtt.removeDir(tmpDir)
    def test_memoryTest_4(self):
        """ --wigglePairs options should be memory clean
        """
        valgrind = mtt.which('valgrind')
        if valgrind is None:
            return
        mtt.makeTempDirParent()
        tmpDir = os.path.abspath(mtt.makeTempDir('memory_4'))
        for maf1, maf2, ref, partner, binlen, valuesDictList in knownValuesWiggles:
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
                    '--wigglePairs=%s:%s' % (ref, partner),
                    '--wiggleBinLength=%d' % binlen,
                    ]
            mtt.recordCommands([cmd], tmpDir)
            mtt.runCommandsS([cmd], tmpDir)
            passed = mtt.noMemoryErrors(os.path.join(tmpDir, 'valgrind.xml'))
            self.assertTrue(passed)
        mtt.removeDir(tmpDir)

class LegitSequencesTest(unittest.TestCase):
    def test_legitSequences(self):
        """ --legitSequences options should produce expected and known output for given inputs
        """
        mtt.makeTempDirParent()
        tmpDir = os.path.abspath(mtt.makeTempDir('legitSequences'))
        i = -1
        for maf1, maf2, totalTrue, totalFalse in knownValues:
            i += 1
            testMaf1 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf1.maf')), 
                                    maf1, g_headers)
            testMaf2 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf2.maf')), 
                                    maf2, g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = []
            cmd += [os.path.abspath(os.path.join(parent, 'test', 'mafComparator')),
                    '--maf1', os.path.abspath(os.path.join(tmpDir, 'maf1.maf')),
                    '--maf2', os.path.abspath(os.path.join(tmpDir, 'maf2.maf')),
                    '--out', os.path.abspath(os.path.join(tmpDir, 'output.xml')),
                    '--samples=1000', '--logLevel=critical',
                    '--legitSequences=%s' % knownValuesLegit[i],
                    ]
            mtt.recordCommands([cmd], tmpDir)
            mtt.runCommandsS([cmd], tmpDir)
            passedTT = totalTrue == getAggregateResult(os.path.abspath(os.path.join(tmpDir, 'output.xml')), 'totalTrue')
            passedTF = totalFalse == getAggregateResult(os.path.abspath(os.path.join(tmpDir, 'output.xml')), 'totalFalse')
            if not (passedTT and passedTF):
                print 'knownValues Test failed on test %d' % i
            self.assertTrue(passedTT and passedTF)
        mtt.removeDir(tmpDir)
    def test_legitSequencesFail(self):
        """ --legitSequences options should fail when sequence length numbers differ from input files
        """
        mtt.makeTempDirParent()
        tmpDir = os.path.abspath(mtt.makeTempDir('legitSequencesFail'))
        i = -1
        for maf1, maf2, totalTrue, totalFalse in knownValues:
            i += 1
            testMaf1 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf1.maf')), 
                                    maf1, g_headers)
            testMaf2 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf2.maf')), 
                                    maf2, g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = []
            cmd += [os.path.abspath(os.path.join(parent, 'test', 'mafComparator')),
                    '--maf1', os.path.abspath(os.path.join(tmpDir, 'maf1.maf')),
                    '--maf2', os.path.abspath(os.path.join(tmpDir, 'maf2.maf')),
                    '--out', os.path.abspath(os.path.join(tmpDir, 'output.xml')),
                    '--samples=1000', '--logLevel=critical',
                    '--legitSequences=%s' % knownValuesLegitFail[i],
                    ]
            mtt.recordCommands([cmd], tmpDir)
            passed = False
            try:
                # we supress the normal stderr output here for the sake of sparing users confusion
                mtt.runCommandsS([cmd], tmpDir, errPipes=[subprocess.PIPE])
            except RuntimeError:
                passed = True
            self.assertTrue(passed)
        mtt.removeDir(tmpDir)
    def test_memoryTest_5(self):
        valgrind = mtt.which('valgrind')
        if valgrind is None:
            return
        mtt.makeTempDirParent()
        tmpDir = os.path.abspath(mtt.makeTempDir('memory_5'))
        i = -1
        for maf1, maf2, totalTrue, totalFalse in knownValues:
            i += 1
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
                    '--legitSequences=%s' % knownValuesLegit[i],
                    ]
            mtt.recordCommands([cmd], tmpDir)
            mtt.runCommandsS([cmd], tmpDir)
            passed = mtt.noMemoryErrors(os.path.join(tmpDir, 'valgrind.xml'))
            self.assertTrue(passed)
        mtt.removeDir(tmpDir)

class NumberOfPairsTest(unittest.TestCase):
    def test_numberOfPairs(self):
        """ --numberOfPairs options should produce expected and known output for given inputs
        """
        mtt.makeTempDirParent()
        tmpDir = os.path.abspath(mtt.makeTempDir('numberOfPairs'))
        i = -1
        for maf1, maf2, totalTrue, totalFalse in knownValues:
            i += 1
            testMaf1 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf1.maf')), 
                                    maf1, g_headers)
            testMaf2 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf2.maf')), 
                                    maf2, g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = []
            cmd += [os.path.abspath(os.path.join(parent, 'test', 'mafComparator')),
                    '--maf1', os.path.abspath(os.path.join(tmpDir, 'maf1.maf')),
                    '--maf2', os.path.abspath(os.path.join(tmpDir, 'maf2.maf')),
                    '--out', os.path.abspath(os.path.join(tmpDir, 'output.xml')),
                    '--samples=1000', '--logLevel=critical',
                    '--numberOfPairs=%s' % knownValuesNumberOfPairs[i],
                    ]
            mtt.recordCommands([cmd], tmpDir)
            mtt.runCommandsS([cmd], tmpDir)
            passedTT = totalTrue == getAggregateResult(os.path.abspath(os.path.join(tmpDir, 'output.xml')), 'totalTrue')
            passedTF = totalFalse == getAggregateResult(os.path.abspath(os.path.join(tmpDir, 'output.xml')), 'totalFalse')
            if not (passedTT and passedTF):
                print 'knownValues Test failed on test %d' % i
            self.assertTrue(passedTT and passedTF)
        mtt.removeDir(tmpDir)
    def test_numberOfPairsFail(self):
        """ --numberOfPairs options should fail for bad input values
        """
        mtt.makeTempDirParent()
        tmpDir = os.path.abspath(mtt.makeTempDir('numberOfPairsFail'))
        i = -1
        for maf1, maf2, totalTrue, totalFalse in knownValues:
            i += 1
            testMaf1 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf1.maf')), 
                                    maf1, g_headers)
            testMaf2 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf2.maf')), 
                                    maf2, g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = []
            cmd += [os.path.abspath(os.path.join(parent, 'test', 'mafComparator')),
                    '--maf1', os.path.abspath(os.path.join(tmpDir, 'maf1.maf')),
                    '--maf2', os.path.abspath(os.path.join(tmpDir, 'maf2.maf')),
                    '--out', os.path.abspath(os.path.join(tmpDir, 'output.xml')),
                    '--samples=1000', '--logLevel=critical',
                    '--numberOfPairs=%s' % knownValuesNumberOfPairsFail[i],
                    ]
            mtt.recordCommands([cmd], tmpDir)
            passed = False
            try:
                mtt.runCommandsS([cmd], tmpDir, errPipes=[subprocess.PIPE])
            except RuntimeError:
                passed = True
            self.assertTrue(passed)
        mtt.removeDir(tmpDir)
    def test_memoryTest_6(self):
        valgrind = mtt.which('valgrind')
        if valgrind is None:
            return
        mtt.makeTempDirParent()
        tmpDir = os.path.abspath(mtt.makeTempDir('memory_6'))
        i = -1
        for maf1, maf2, totalTrue, totalFalse in knownValues:
            i += 1
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
                    '--numberOfPairs=%s' % knownValuesNumberOfPairs[i],
                    ]
            mtt.recordCommands([cmd], tmpDir)
            mtt.runCommandsS([cmd], tmpDir)
            passed = mtt.noMemoryErrors(os.path.join(tmpDir, 'valgrind.xml'))
            self.assertTrue(passed)
        mtt.removeDir(tmpDir)

class KnownValuesTest(unittest.TestCase):
    # knownValues contains quad-tuples,
    # maf1, maf2, totalTrue (comparing maf1 as fileA to maf2), 
    # totalFalse (comparing maf1 as fileA to maf2)
    def test_knownValues(self):
        """ mafComparator should return correct results for hand-calculable problems
        """
        mtt.makeTempDirParent()
        tmpDir = os.path.abspath(mtt.makeTempDir('knownValues'))
        i = -1
        for maf1, maf2, totalTrue, totalFalse in knownValues:
            i += 1
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
            if not (passedTT and passedTF):
                print 'knownValues Test failed on test %d' % i
            self.assertTrue(passedTT and passedTF)
        mtt.removeDir(tmpDir)
    def test_memory_2(self):
        """ mafComparator should be memory clean for known values
        """
        valgrind = mtt.which('valgrind')
        if valgrind is None:
            return
        mtt.makeTempDirParent()        
        tmpDir = os.path.abspath(mtt.makeTempDir('memory_2'))
        for maf1, maf2, totalTrue, totalFalse in knownValues:
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
    def test_bedParsing(self):
        """ mafComparator should parse a bed file and use the intervals for testing
        """
        mtt.makeTempDirParent()
        tmpDir = os.path.abspath(mtt.makeTempDir('bedParsing'))
        for maf1, maf2, bed, totalTrue, totalTrueInInterval in knownValuesBed:
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
        for maf1, maf2, bed, totalTrue, totalTrueInInterval in knownValuesBed:
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
    def test_seedTesting(self):
        """ mafComparator should have replicable runs via the --seed command
        """
        mtt.makeTempDirParent()
        tmpDir = os.path.abspath(mtt.makeTempDir('seedTesting'))
        parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        for maf1, maf2  in knownValuesSeed:
            testMaf1 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf1.maf')), 
                                    maf1, g_headers)
            testMaf2 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf2.maf')), 
                                    maf2, g_headers)
            cmd = [os.path.abspath(os.path.join(parent, 'test', 'mafComparator')),
                   '--mafFile1', os.path.abspath(os.path.join(tmpDir, 'maf1.maf')), 
                   '--mafFile2', os.path.abspath(os.path.join(tmpDir, 'maf2.maf')),
                   '--outputFile', os.path.join(tmpDir, 'output.xml'),
                   '--samples=10', '--logLevel=critical']
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
        parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        for maf1, maf2  in knownValuesSeed:
            testMaf1 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf1.maf')), 
                                    maf1, g_headers)
            testMaf2 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf2.maf')), 
                                    maf2, g_headers)
            cmd = mtt.genericValgrind(tmpDir)
            cmd += [os.path.abspath(os.path.join(parent, 'test', 'mafComparator')),
                    '--maf1', os.path.abspath(os.path.join(tmpDir, 'maf1.maf')), 
                    '--maf2', os.path.abspath(os.path.join(tmpDir, 'maf2.maf')),
                    '--out', os.path.join(tmpDir, 'output.xml'),
                    '--samples=10', '--logLevel=critical']
            mtt.recordCommands([cmd], tmpDir)
            mtt.runCommandsS([cmd], tmpDir)
            self.assertTrue(mtt.noMemoryErrors(os.path.join(tmpDir, 'valgrind.xml')))
            for i in xrange(0, 4):
                mtt.recordCommands([cmd], tmpDir)
                mtt.runCommandsS([cmd], tmpDir)
                self.assertTrue(mtt.noMemoryErrors(os.path.join(tmpDir, 'valgrind.xml')))
        mtt.removeDir(tmpDir)
class NearTests(unittest.TestCase):
    def test_nearSimple(self):
        """ mafComparator should return correct results for hand-calculable problems that use the --near=0 option
        """
        mtt.makeTempDirParent()
        tmpDir = os.path.abspath(mtt.makeTempDir('nearSimple'))
        i = 0
        for maf1, maf2, totalTrue, totalFalse in knownValues:
            i += 1
            testMaf1 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf1.maf')), 
                                    maf1, g_headers)
            testMaf2 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf2.maf')), 
                                    maf2, g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = [os.path.abspath(os.path.join(parent, 'test', 'mafComparator')),
                   '--maf1', os.path.abspath(os.path.join(tmpDir, 'maf1.maf')),
                   '--maf2', os.path.abspath(os.path.join(tmpDir, 'maf2.maf')),
                   '--out', os.path.abspath(os.path.join(tmpDir, 'output.xml')),
                   '--near=0', 
                   '--samples=1000', '--logLevel=critical',
                   ]
            mtt.recordCommands([cmd], tmpDir)
            mtt.runCommandsS([cmd], tmpDir)
            passedTT = totalTrue == getAggregateResult(os.path.abspath(os.path.join(tmpDir, 'output.xml')), 'totalTrue')
            passedTF = totalFalse == getAggregateResult(os.path.abspath(os.path.join(tmpDir, 'output.xml')), 'totalFalse')
            self.assertTrue(passedTT and passedTF)
        mtt.removeDir(tmpDir)
    def test_near(self):
        """ mafComparator should return correct results for hand-calculable problems that use the --near option
        """
        mtt.makeTempDirParent()
        tmpDir = os.path.abspath(mtt.makeTempDir('near'))
        i = 0
        for maf1, maf2, near, totalTrue, totalFalse in knownValuesNear:
            i += 1
            testMaf1 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf1.maf')), 
                                    maf1, g_headers)
            testMaf2 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf2.maf')), 
                                    maf2, g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = [os.path.abspath(os.path.join(parent, 'test', 'mafComparator')),
                   '--maf1', os.path.abspath(os.path.join(tmpDir, 'maf1.maf')),
                   '--maf2', os.path.abspath(os.path.join(tmpDir, 'maf2.maf')),
                   '--out', os.path.abspath(os.path.join(tmpDir, 'output.xml')),
                   '--near=%d' % near, 
                   '--samples=1000', '--logLevel=critical',
                   ]
            mtt.recordCommands([cmd], tmpDir)
            mtt.runCommandsS([cmd], tmpDir)
            passedTT = totalTrue == getAggregateResult(os.path.abspath(os.path.join(tmpDir, 'output.xml')), 'totalTrue')
            passedTF = totalFalse == getAggregateResult(os.path.abspath(os.path.join(tmpDir, 'output.xml')), 'totalFalse')
            self.assertTrue(passedTT and passedTF)
        mtt.removeDir(tmpDir)
    def test_memory_3(self):
        """ mafComparator should be memory clean for known values using --near option
        """
        valgrind = mtt.which('valgrind')
        if valgrind is None:
            return
        mtt.makeTempDirParent()        
        tmpDir = os.path.abspath(mtt.makeTempDir('memory_3'))
        for maf1, maf2, near, totalTrue, totalFalse in knownValuesNear:
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
                    '--near=%d' % near, 
                    '--samples=1000', '--logLevel=critical',
                    ]
            mtt.recordCommands([cmd], tmpDir)
            mtt.runCommandsS([cmd], tmpDir)
            passed = mtt.noMemoryErrors(os.path.join(tmpDir, 'valgrind.xml'))
            self.assertTrue(passed)
        mtt.removeDir(tmpDir)
        
class CuTestTests(unittest.TestCase):
    ##################################################
    # DISABLED DUE to excessive amout of time this takes to run
    ##################################################
    def gtest_CuTestTests(self):
        """ Yo dawg, I heard you liked unit tests so I put some unit tests in your unit test so now you can unit test when you unit test.
        """
        mtt.makeTempDirParent()
        tmpDir = os.path.abspath(mtt.makeTempDir('CuTestAllTests'))
        cmd = mtt.genericValgrind(tmpDir)
        cmd += [os.path.abspath(os.path.join(os.curdir, 'test', 'allTests'))]
        outpipes = [os.path.join('/dev', 'null')]
        mtt.recordCommands([cmd], tmpDir)
        mtt.runCommandsS([cmd], tmpDir, outPipes=outpipes)
        self.assertTrue(mtt.noMemoryErrors(os.path.join(tmpDir, 'valgrind.xml')))
        mtt.removeDir(tmpDir)

if __name__ == '__main__':
    unittest.main()
