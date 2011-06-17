import unittest
import os, sys
myBinDir = os.path.normpath(os.path.dirname(sys.argv[0]))
sys.path.append(myBinDir + "/../../..")
os.environ["PATH"] = myBinDir + "/../../../../bin:" + os.environ["PATH"]
from sonLib.bioio import logger

class VerifyMafJoinOutput( unittest.TestCase ):
    import os
    if not os.path.exists('temp_testFiles'):
        os.mkdir('temp_testFiles')
    def test_knownValues ( self ):
        """mafJoin output should produce no missing when run in eval_MAFComparator"""
        import subprocess
        import sys
        childA = 'simHuman'
        childB = 'simGorilla'
        parent = 'sG-sH-sC'
        cmd = ['mafJoin', '-treelessRoot1=%s' % parent, '-treelessRoot2=%s' % parent,
               parent, '-maxBlkWidth=1000',
               '-multiParentDropped=temp_testFiles/%s.dropped.tab' % parent,
               'evolverMAFs/%s.%s.avg30.maf' % ( childA, parent ),
               'evolverMAFs/%s.%s.avg30.maf' % ( childB, parent ),
               'temp_testFiles/%s.mafJoin.maf' %parent ]
        logger.info("run: " + " ".join(cmd))
        p = subprocess.Popen( cmd )
        p.wait()

        dropped = open('temp_testFiles/%s.dropped.tab' % parent).read()
        self.assertEqual('maf\tdroppedSeq\tdroppedStart\tdroppedEnd\n', dropped )
        
        cmd = ['eval_MAFComparator', '--mAFFile1=evolverMAFs/%s.%s.avg30.maf' % ( childA, parent ),
               '--mAFFile2=temp_testFiles/%s.mafJoin.maf' % parent,
               '--outputFile=temp_testFiles/%s.%s.compare.xml' % ( childA, parent ),
               '--sampleNumber=100000000', '--ultraVerbose']
        missingFile = open('temp_testFiles/%s.%s.missing.tab' % ( childA, parent ),'w')
        logger.info("run: " + " ".join(cmd))
        p = subprocess.Popen( cmd, stderr=missingFile )
        p.wait()

        missing = open( 'temp_testFiles/%s.%s.missing.tab' % ( childA, parent )).read()
        emptyMissing  = '# Comparing evolverMAFs/%s.%s.avg30.maf to temp_testFiles/%s.mafJoin.maf\n' % ( childA, parent, parent )
        emptyMissing += '# seq1\tabsPos1\torigPos1\tseq2\tabsPos2\torigPos2\n'
        line1 = len( emptyMissing )
        emptyMissing += '# Comparing temp_testFiles/%s.mafJoin.maf to evolverMAFs/%s.%s.avg30.maf\n' % ( parent, childA, parent )
        emptyMissing += '# seq1\tabsPos1\torigPos1\tseq2\tabsPos2\torigPos2\n'
        bothLines = len( emptyMissing )

        self.assertEqual( emptyMissing[ 0:line1 ], missing[ 0:line1 ] ) # this should be true no matter what, it's just the two headers.
        self.assertEqual( emptyMissing[ 0:bothLines ], missing[ 0:bothLines ] )
        # it's fine if there's more to the file after this point,
        # what matters is if there are any missing items in the first set of comparisons.

# FIXME: should be controllable from the command line
import logging
logger.setLevel(logging.INFO)
if __name__ == "__main__":
    unittest.main()
