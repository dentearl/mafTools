import unittest
import os, sys
myBinDir = os.path.normpath(os.path.dirname(sys.argv[0]))
sys.path.append(myBinDir + "/../../..")
os.environ["PATH"] = myBinDir + "/../../../../bin:" + os.environ["PATH"]
from sonLib.bioio import logger

class VerifyMafJoinOutput( unittest.TestCase ):
    import os
    if not os.path.exists('temp_offsetTests'):
        os.mkdir('temp_offsetTests')
    def test_knownValues ( self ):
        """mafJoin output should produce no missing when run in eval_MAFComparator"""
        import subprocess
        import sys
        cmd = ['mafJoin', '-treelessRoot1=sMouse-sRat', '-treelessRoot2=sMouse-sRat',
               'sMouse-sRat', '-maxBlkWidth=1000',
               '-multiParentDropped=temp_offsetTests/sMouse-sRat.dropped.tab',
               'evolverMAFs/simRat.sMouse-sRat.maf',
               'evolverMAFs/simMouse.sMouse-sRat.maf',
               'temp_offsetTests/sMouse-sRat.mafJoin.maf']
        logger.info("run: " + " ".join(cmd))
        p = subprocess.Popen( cmd )
        p.wait()
        dropped = open('temp_offsetTests/sMouse-sRat.dropped.tab').read()

        self.assertEqual('maf\tdroppedSeq\tdroppedStart\tdroppedEnd\n', dropped )
        
        cmd = ['eval_MAFComparator', '--mAFFile1=evolverMAFs/simRat.sMouse-sRat.maf',
               '--mAFFile2=temp_offsetTests/sMouse-sRat.mafJoin.maf',
               '--outputFile=temp_offsetTests/simRat.sMouse-sRat.compare.xml',
               '--sampleNumber=100000000', '--ultraVerbose']
        missingFile = open('temp_offsetTests/simRat.sMouse-sRat.missing.tab','w')
        logger.info("run: " + " ".join(cmd))
        p = subprocess.Popen( cmd, stderr=missingFile )
        p.wait()

        missing = open( 'temp_offsetTests/simRat.sMouse-sRat.missing.tab' ).read()
        emptyMissing  = '# Comparing evolverMAFs/simRat.sMouse-sRat.maf to temp_offsetTests/sMouse-sRat.mafJoin.maf\n'
        emptyMissing += '# seq1\tpos1\tseq2\tpos2\n'
        emptyMissing += '# Comparing temp_offsetTests/sMouse-sRat.mafJoin.maf to evolverMAFs/simRat.sMouse-sRat.maf\n'
        emptyMissing += '# seq1\tpos1\tseq2\tpos2\n'

        self.assertEqual( emptyMissing[0:111], missing[0:111] ) # this should be true no matter what, it's just the two headers.
        self.assertEqual( emptyMissing[0:225], missing[0:225] )

# FIXME: should be controllable from the command line
import logging
logger.setLevel(logging.INFO)
if __name__ == "__main__":
    unittest.main()
