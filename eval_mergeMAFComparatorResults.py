#!/usr/bin/env python

"""Script for computing alignments for a reconstruction problem.
"""

import xml.etree.ElementTree as ET
import os

from sonLib.bioio import getBasicOptionParser
from sonLib.bioio import parseBasicOptions
from sonLib.bioio import logger

def main():
    ##########################################
    #Construct the arguments.
    ##########################################    
    
    parser = getBasicOptionParser("usage: %prog [options]", "%prog 0.1")
    
    parser.add_option("--results1", dest="results1", 
                      help="File containing the first XML formatted MAF comparsion results.")
    
    parser.add_option("--results2", dest="results2", 
                      help="File containing the second XML formatted MAF comparsion results.")
    
    parser.add_option("--outputFile", dest="outputFile", type="string",
                      help="The file to put the aggregated results in (can be the same as either of the two inputs)")
        
    options, args = parseBasicOptions(parser)
        
    assert len(args) == 0
    logger.info("Parsed arguments")
    
    assert options.results1 != None
    assert options.results2 != None
    assert options.outputFile != None
    
    ##########################################
    #Do the merging.
    ##########################################
    
    resultsTree1 = ET.parse(options.results1).getroot()
    resultsTree2 = ET.parse(options.results2).getroot()
    
    homologyTestsList1 = resultsTree1.findall("homology_tests")
    homologyTestsList2 = resultsTree2.findall("homology_tests")
    
    assert len(homologyTestsList1) == len(homologyTestsList2)
    for i in xrange(len(homologyTestsList1)):
        homologyTests1 = homologyTestsList1[i]
        homologyTests2 = homologyTestsList2[i]
        #Do merge for tests in both sets
        for homologyTest1 in homologyTests1.findall("homology_test"):
            for homologyTest2 in homologyTests2.findall("homology_test"):
                if homologyTest1.attrib["sequenceA"] == homologyTest2.attrib["sequenceA"] and \
                homologyTest1.attrib["sequenceB"] == homologyTest2.attrib["sequenceB"]:
                    
                    homologyTest1.attrib["totalTests"] = str(int(homologyTest1.attrib["totalTests"]) + int(homologyTest2.attrib["totalTests"]))
                    homologyTest1.attrib["totalTrue"] = str(int(homologyTest1.attrib["totalTrue"]) + int(homologyTest2.attrib["totalTrue"]))
                    homologyTest1.attrib["totalFalse"] = str(int(homologyTest1.attrib["totalFalse"]) + int(homologyTest2.attrib["totalFalse"]))
                    if int(homologyTest1.attrib["totalTests"]) != 0:
                        homologyTest1.attrib["average"] = str(float(homologyTest1.attrib["totalTrue"]) /  float(homologyTest1.attrib["totalTests"]))
                    
                    for singleHomologyTest2 in homologyTest2.findall("single_homology_test"):
                        homologyTest1.insert(0, singleHomologyTest2)
        
        #Now add in tests not in the intersection of the results
        l = []
        for homologyTest2 in homologyTests2.findall("homology_test"):
            for homologyTest1 in homologyTests1.findall("homology_test"):
                if homologyTest1.attrib["sequenceA"] == homologyTest2.attrib["sequenceA"] and \
                    homologyTest1.attrib["sequenceB"] == homologyTest2.attrib["sequenceB"]:
                    break
            else:
                l.append(homologyTest2)
        
        for homologyTest2 in l:
            homologyTests1.insert(0, homologyTest2)
            
        #Now recalculate the totals
        totalTrue = 0
        totalFalse = 0
        totalTests = 0
        for homologyTest1 in homologyTests1.findall("homology_test"):
            totalTests += int(homologyTest1.attrib["totalTests"])
            totalTrue += int(homologyTest1.attrib["totalTrue"])
            totalFalse += int(homologyTest1.attrib["totalFalse"])
        assert totalTrue + totalFalse == totalTests
        homologyTests1.attrib["totalTests"] = str(totalTests)
        homologyTests1.attrib["totalTrue"] = str(totalTrue)
        homologyTests1.attrib["totalFalse"] = str(totalFalse)
        homologyTests1.attrib["average"] = str(float(totalTrue)/(totalTests + 0.01))
        
    
    #Write to the results file.
    fileHandle = open(options.outputFile, 'w')
    tree = ET.ElementTree(resultsTree1)
    tree.write(fileHandle)
    fileHandle.close()
    

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    _test()
    main()
