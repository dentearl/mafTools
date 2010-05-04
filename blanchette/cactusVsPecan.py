import os
from sonLib.bioio import getTempDirectory
from sonLib.bioio import logger
from sonLib.bioio import system
from sonLib.bioio import TestStatus

"""
For each sub-directory containing a simulation:
Make a temporary directory to do the work in
Convert the 'true' alignment MFA file to a MAF file.
Run cactus, generate a netDisk, convert this to a MAF file.
Run Pecan, generate a MFA file, convert this to a MAF file.
Run eval_MAFComparator on the true alignment's MAF and the Cactus MAF
Run eval_MAFComparator on the true alignment's MAF and the Pecan MAF
For each of the two results files,  if there already exists some result, merge it with the new results using eval_mergeMAFComparatorResults
else copy the results to the final output files.
Empty the temp dir for next time.
"""

blanchettePath = os.path.join(TestStatus.getPathToDataSets(), "blanchettesSimulation")
newickTreeString = "((((HUMAN:0.006969, CHIMP:0.009727):0.025291, BABOON:0.044568):0.11,(MOUSE:0.072818, RAT:0.081244):0.260342):0.023260,((DOG:0.07, CAT:0.07):0.087381,(PIG:0.06, COW:0.06):0.104728):0.04);"
species = ("HUMAN", "CHIMP", "BABOON", "MOUSE", "RAT", "DOG", "CAT", "PIG", "COW")
tempDir = getTempDirectory(os.getcwd())
simulationNumber = 1
sampleNumber = 1000000
cactusAggregateResults = "cactus.xml"
pecanAggregateResults = "pecan.xml"

for simulation in xrange(simulationNumber):
    logger.info("Starting simuation number %i" % simulation)
    simulationPath = os.path.join(blanchettePath, "%.2i.job" % simulation)
    sequences = [ os.path.join(simulationPath, speci) for speci in species ] #Joke
    trueMFA = os.path.join(simulationPath, "true.mfa")
    trueMAF = os.path.join(tempDir, "true.maf")
    logger.info("Got the sim path %s, sequence files: %s, the true MFA: %s and the true MAF %s" % (simulationPath, " ".join(sequences), trueMFA, trueMAF))
    system("eval_MFAToMAF --mFAFile %s --outputFile %s --logLevel DEBUG" % (trueMFA, trueMAF))
    logger.info("Converted the true MFA to MAF")
    jobTree = os.path.join(tempDir, "jobTree")
    netDisk = os.path.join(tempDir, "netDisk")
    system("jobTree.py --logLevel DEBUG --jobTree %s --command \"cactus_workflow.py --job JOB_FILE --speciesTree '%s' --netDisk %s %s --setupAndBuildAlignments\"" % (jobTree, newickTreeString, netDisk, " ".join(sequences)))
    system("jobTreeStatus.py --jobTree %s --failIfNotComplete --logLevel DEBUG" % (jobTree,))
    logger.info("Ran cactus tree apparently okay")
    cactusMAF = os.path.join(tempDir, "cactus.maf")
    system("cactus_MAFGenerator --logLevel DEBUG --netDisk %s --outputFile %s --netName 0" % (netDisk, cactusMAF))
    logger.info("Got the MAFS from Cactus tree.")
    #pecanMFA = os.path.join(tempDir, "pecan.mfa")
    #system("java bp.pecan.Pecan -A -E '%s' -F %s -G %s" % (newickTreeString, " ".join(sequences), pecanMFA))
    #logger.info("Ran Pecan okay")
    #pecanMAF = os.path.join(tempDir, "pecan.maf")
    #system("eval_MFAToMAF --logLevel DEBUG --mFAFile %s --outputFile %s " % (pecanMFA, pecanMAF))
    #logger.info("Got the MAF from the Pecan MFA")
    cactusResults = os.path.join(tempDir, "cactus.xml")
    system("eval_MAFComparator --logLevel DEBUG --mAFFile1 %s --mAFFile2 %s --outputFile %s --sampleNumber %i" % (trueMAF, cactusMAF, cactusResults, sampleNumber))
    logger.info("Ran the true MAF to cactus MAF comparison")
    #pecanResults = os.path.join(tempDir, "pecan.xml")
    #system("eval_MAFComparator --logLevel DEBUG --mAFFile1 %s --mAFFile2 %s --outputFile %s --sampleNumber %i" % (trueMAF, pecanMAF, pecanResults, sampleNumber))
    #logger.info("Ran the true MAF to pecan MAF comparison")
    if simulation > 0:
        system("eval_mergeMAFComparatorResults.py --logLevel DEBUG --results1 %s --results2 %s --outputFile %s" % (cactusAggregateResults, cactusResults, cactusAggregateResults))
        #system("eval_mergeMAFComparatorResults.py --logLevel DEBUG --results1 %s --results2 %s --outputFile %s" % (pecanAggregateResults, pecanResults, pecanAggregateResults))
        logger.info("Merged together the results of previous simulations with these ones")
    else:
        system("mv %s %s" % (cactusResults, cactusAggregateResults))
        #system("mv %s %s" % (pecanResults, pecanAggregateResults))
    
    #Run cactus tree stats
    cactusStats = os.path.join("./", "cactus%s.xml" % simulation)
    system("cactus_treeStats --logLevel DEBUG --netDisk %s --outputFile %s --netName 0" % (netDisk, cactusStats))  
    
    system("rm -rf %s/*" % tempDir) 
    logger.info("Cleaned up the temp dir for the next run")
   
system("rm -rf %s" % tempDir)
logger.info("Cleaned up and am done")