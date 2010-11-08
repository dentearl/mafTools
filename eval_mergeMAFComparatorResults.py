#!/usr/bin/env python

"""Script for computing alignments for a reconstruction problem.
"""

import xml.etree.ElementTree as ET
import os
import sys

from xml.dom import minidom

from sonLib.bioio import getBasicOptionParser
from sonLib.bioio import parseBasicOptions
from sonLib.bioio import logger

# Taken from:
# http://blog.doughellmann.com/2010/03/pymotw-creating-xml-documents-with.html
def prettify(elem):
    """Return a pretty-printed XML string for the Element.
    """
    rough_string = ET.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="\t")

def main():
    ##########################################
    #Construct the arguments.
    ##########################################    
    
    parser = getBasicOptionParser("usage: %prog [options]", "%prog 0.1")
    
    parser.add_option("--results1", dest="results1", 
                      help="File containing the first XML formatted MAF comparsion results.")
    
    parser.add_option("--results2", dest="results2", 
                      help="File containing the second XML formatted MAF comparsion results.")

    parser.add_option("--inputFile", dest="inputFile",
                      help="File containing a list of XML formatted MAF comparison results.")

    parser.add_option("--inputDir", dest="inputDir",
                      help="Directory containing the XML formatted MAF comparison results.")
    
    parser.add_option("--outputFile", dest="outputFile", type="string",
                      help="The file to put the aggregated results in (can be the same as either of the two inputs)")
        
    options, args = parseBasicOptions(parser)
        
    assert len(args) == 0
    logger.info("Parsed arguments")
    
    assert (options.results1 != None and options.results2 != None) or (options.inputDir != None) or (options.inputFile != None)
    assert options.outputFile != None
    
    ##########################################
    #Do the merging.
    ##########################################

    xmlList = []
    if options.inputDir != None:
        try:
            fileList = os.listdir(options.inputDir)
        except:
            print >> sys.stderr, "Error: Can't open dir '%s'" % options.inputDir
            sys.exit(-1)

        for file in fileList:
            if file.endswith(".xml"):
                if options.inputDir.startswith("/"):
                    xmlList.append(options.inputDir + "/" + file)
                else:
                    xmlList.append(os.getcwd() + "/" + options.inputDir + "/" + file)
    elif options.inputFile != None:
        try:
            f = open(options.inputFile)
        except:
            print >> sys.stderr, "Error: Can't open file '%s'" % options.inputFile
            sys.exit(-1)
        for line in f:
            line = line.rstrip()
            xmlList.append(line)
        f.close()
    elif options.results1 != None and options.results2 != None:
        xmlList.append(options.results1)
        xmlList.append(options.results2)
    else:
        print >> sys.stderr, "Error: Need to specify input"
        sys.exit(-1)

    baseFile = xmlList.pop()
    resultsTree1 = ET.parse(baseFile).getroot()
    homologyTestsList1 = resultsTree1.findall("homology_tests")
    for i in xrange(len(homologyTestsList1)):
        homologyTests1 = homologyTestsList1[i]
        for homologyTest1 in homologyTests1.findall("homology_test"):
            for singleHomologyTest1 in homologyTest1.findall("single_homology_test"):
                singleHomologyTest1.attrib["srcFile"] = str(baseFile)

    for xmlFile in xmlList:
        resultsTree2 = ET.parse(xmlFile).getroot()
        homologyTestsList2 = resultsTree2.findall("homology_tests")
    
        assert len(homologyTestsList1) == len(homologyTestsList2)
        for i in xrange(len(homologyTestsList1)):
            homologyTests1 = homologyTestsList1[i]
            homologyTests2 = homologyTestsList2[i]
            assert(homologyTests1.attrib["near"] == homologyTests2.attrib["near"])
            
            def mergeAggregateResults(tag1, tag2):
                def sum(i, j):
                    i.attrib["totalTests"] = str(int(i.attrib["totalTests"]) + int(j.attrib["totalTests"]))
                    i.attrib["totalTrue"] = str(int(i.attrib["totalTrue"]) + int(j.attrib["totalTrue"]))
                    i.attrib["totalFalse"] = str(int(i.attrib["totalFalse"]) + int(j.attrib["totalFalse"]))
                    if int(i.attrib["totalTests"]) != 0:
                        i.attrib["average"] = str(float(i.attrib["totalTrue"]) / float(i.attrib["totalTests"]))
                assert tag1.tag == "aggregate_results"
                assert tag2.tag == "aggregate_results"
                sum(tag1.find("all"), tag2.find("all"))
                sum(tag1.find("both"), tag2.find("both"))
                sum(tag1.find("A"), tag2.find("A"))
                sum(tag1.find("B"), tag2.find("B"))
                sum(tag1.find("neither"), tag2.find("neither"))
            
            #Do merge for tests in both sets
            for homologyTest1 in homologyTests1.find("homology_pair_tests").findall("homology_test"):
                for homologyTest2 in homologyTests2.find("homology_pair_tests").findall("homology_test"):
                    if homologyTest1.attrib["sequenceA"] == homologyTest2.attrib["sequenceA"] and homologyTest1.attrib["sequenceB"] == homologyTest2.attrib["sequenceB"]:
                        mergeAggregateResults(homologyTest1.find("aggregate_results"), homologyTest2.find("aggregate_results"))
                        for singleHomologyTest2 in homologyTest2.find("single_homology_tests").findall("single_homology_test"):
                            homologyTest1.find("single_homology_tests").insert(0, singleHomologyTest2)
        
            #Now add in tests not in the intersection of the results
            l = []
            for homologyTest2 in homologyTests2.find("homology_pair_tests").findall("homology_test"):
                for homologyTest1 in homologyTests1.find("homology_pair_tests").findall("homology_test"):
                    if homologyTest1.attrib["sequenceA"] == homologyTest2.attrib["sequenceA"] and homologyTest1.attrib["sequenceB"] == homologyTest2.attrib["sequenceB"]:
                        break
                else:
                    l.append(homologyTest2)
        
            for homologyTest2 in l:
                homologyTests1.find("homology_pair_tests").insert(0, homologyTest2)
            
            #Now recalculate the totals
            mergeAggregateResults(homologyTests1.find("aggregate_results"), homologyTests2.find("aggregate_results"))
        
    #Write to the results file.
    fileHandle = open(options.outputFile, 'w')
#    tree = ET.ElementTree(resultsTree1)
#    tree.write(fileHandle)
    fileHandle.write(prettify(resultsTree1))
    fileHandle.close()

    return
    
def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    _test()
    main()
