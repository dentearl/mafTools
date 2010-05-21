"""Test for eval_mAFJoiner. 

Let A be a MAF file containing a tree for each MAF block.
We define the following properties for A.

Let a be a MAF block
in A. Let T(a) be the 'block tree' for a. T(a) is composed of 'leaf nodes', with no descendants and
'internal nodes', which have 1 or more descendants. The root node of T(a) is the unique most
ancestral node of the tree. For each node in T(a) there is a corresponding label in the tree and row in the MAF block. 
For each node n in T(a) except the root node, there is an associated ancestral branch, this branch must 
be labelled with an associated branch length. 

Each node in T(a) also has an associated 'species' label. Let T'(A) be the 'species tree' associated with the genomes in A. 
We call the nodes in T'(A) 'species nodes', reusing the same terminology as with the MAF block trees.
The tree T(a) for any a in A may 'disagree' with T'(A), due to duplication, loss and incomplete lineage sorting, however
if x is the root species nodes of T'(A) then for any node n with species label x in any tree T(a), n 
must be a root node, i.e. no coalescence can occur that is more ancestral than the coalesence of the root species.

A MAF file with the above properties is called 'tree normal'. Again let A be a tree normal MAF file.
Let i represent a column in a MAF block a in A. 
Let F(i) be the set of sequence positions in i (i.e. excluding gaps).
Let T''(i) be the tree for column i, which is an induced subtree of T(a) excluding nodes that are gapped in the alignment. 
Each node in T''(i) therefore has an associated sequence position in F(i).
Let S(A) be the species in T'(A).

Let i and j now represent two columns from two different MAF files. i and j can be 'joined' 
if the intersection of F(i) and F(j) is non-empty to form a larger
alignment column k of sequence positions that are all derived from a common ancestor.
The join to form a column k is 'simple' if T''(k) is a tree composed of T''(i)
and T''(j) merged at a single common node (this needs more explanation). 
A 'complete join' of two tree normal MAF files is a third tree normal file whose columns contain the 'closure' (wrong word) of joins of 
the columns of A and B.

Let A and B be two tree normal MAF files. We say that A and B are 'compatible' iff (1) the intersection of S(A) and S(B) is { x } (so x is the 
unique shared species), and (2) x is the root of either or both of T(A) and T(B).

Theorem 1:  For two compatible tree normal MAF files there exists a (unique) 
complete join of A and B for which all joins are simple.

mafJoiner performs a complete join on two compatible tree normal MAF files.
"""

import os,sys
myBinDir = os.path.normpath(os.path.dirname(sys.argv[0]))
sys.path.append(myBinDir + "/../..")
os.putenv("PATH", myBinDir + "/../../../bin:" + os.getenv("PATH"))

import unittest
from sonLib.bioio import getTempFile 
from sonLib.bioio import logger
from sonLib.bioio import system

class MafJoinTests(unittest.TestCase):

    def getTestTempFile(self, suffix):
        tempDir = "tmp"
        if not os.path.exists(tempDir):
            os.makedirs(tempDir)
        tempFile = tempDir + "/" + self.id() + "." + suffix
        if os.path.exists(tempFile):
            os.unlink(tempFile)
        return tempFile
    
    def writeMAF(self, maf, suffix):
        """Takes one of the above MAF strings and writes it to temp file.
        """
        tempFile = self.getTestTempFile(suffix)
        fileHandle = open(tempFile, 'w')
        fileHandle.write("##maf version=1\n")
        # skip first blank line
        for line in maf.split("\n")[1:]:
            fileHandle.write(line.strip() + "\n")
        fileHandle.close()
        return tempFile

    def runEvalMafJoiner(self, ref, mAFFileA, mAFFileB, outputMAFFile):
        system("mafJoin %s %s %s %s" \
               % (ref, mAFFileA, mAFFileB, outputMAFFile))
        logger.info("Ran eval-maf joiner okay")

    def compareExpectedAndRecieved(self, expected, recieved):
        """Checks two MAFs are equivalent.
        """
        system("diff -u %s %s" % (expected, recieved))

    def mafJoinTest(self, ref, mafA, mafB, mafC):
        """Writes out mafA and mafB to temp files.
        Runs eval_mAFJoiner.
        Parses the output and compares it to mafC.
        """
        tempFileA = self.writeMAF(mafA, "A.maf")
        tempFileB = self.writeMAF(mafB, "B.maf")
        tempFileC = self.writeMAF(mafC, "C.maf")
        tempOutputFile = self.getTestTempFile("out.maf")
        logger.info("Got inputs to test")
        self.runEvalMafJoiner(ref, tempFileA, tempFileB, tempOutputFile)
        self.compareExpectedAndRecieved(tempFileC, tempOutputFile)
        logger.info("Ran test apparently okay")

    def setUp(self):
        unittest.TestCase.setUp(self)
    
    def tearDown(self):
        unittest.TestCase.tearDown(self)
        
    def testJoin1(self):
        """Simple non-dup join. Shows splitting of blocks.
        Note the join process should maintain the ordering of rows in columns.
        """
        A = """
        a score=50.0 tree='(hg18.chr7:0.1,mm4.chr6:0.15)rn3.chr4;'
        s hg18.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
        s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
        s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
        """
        B = """
        a score=5.0 tree='(panTro1.chr6:0.3,baboon:0.2)hg18.chr7;'
        s panTro1.chr6 28741140 33 + 161576975 AAAGGGAATGTTAACCAAATGAATTGTCTCTTA
        s baboon         116834 33 +   4622798 AAAGGGAATGTTAACCAAATGAGTTGTCTCTTA
        s hg18.chr7    27578828 33 + 158545518 AAAGGGAATGTTAACCAAATGAATTGTCTCTTA
        """
        C = """
        a score=0.000000 tree='((panTro1.chr6:0.3,baboon:0.2)hg18.chr7:0.1,mm4.chr6:0.15)rn3.chr4;'
        s panTro1.chr6  28741140 33 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTA-----
        s baboon.baboon   116834 33 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTA-----
        s hg18.chr7     27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
        s mm4.chr6      53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
        s rn3.chr4      81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
        """
        self.mafJoinTest("hg18", A, B, C)
    
    def testJoin2(self):
        """Simple non-dup join. Shows ordering of inserts, first from A then from B.
        """
        A = """
        a score=10.0 tree='(hg18.chr7:0.1,mm4.chr6:0.15)rn3.chr4;'
        s hg18.chr7    27699739 3 + 158545518 T---GA
        s mm4.chr6     53303881 6 + 151104725 TAAAGA
        s rn3.chr4     81444246 6 + 187371129 taagga
        """
        B = """
        a score=1000.0 tree='(panTro1.chr6:0.3,baboon:0.2)hg18.chr7;'
        s panTro1.chr6 28862317 6 + 161576975 TAAAGA
        s baboon         241163 3 +   4622798 T---GA
        s hg18.chr7    27699739 3 + 158545518 T---GA
        """
        C = """
        a score=0.000000 tree='((panTro1.chr6:0.3,baboon:0.2)hg18.chr7:0.1,mm4.chr6:0.15)rn3.chr4;'
        s panTro1.chr6  28862317 6 + 161576975 TAAA---GA
        s baboon.baboon   241163 3 +   4622798 T------GA
        s hg18.chr7     27699739 3 + 158545518 T------GA
        s mm4.chr6      53303881 6 + 151104725 T---AAAGA
        s rn3.chr4      81444246 6 + 187371129 t---aagga
        """
        self.mafJoinTest("hg18", A, B, C)
    
    def testJoin3(self):
        """Simple non-dup join. Shows merging of blocks. Shows merging of trees at a 
        birfurcation.
        """
        A = """
        a score=5.0 tree='(mm4.chr6:0.15)hg18.chr7;'
        s hg18.chr7    27707221 2 + 158545518 gc
        s mm4.chr6     53310102 2 - 151104725 AC
        
        a score=5.0 tree='(mm4.chr6:0.15)hg18.chr7;'
        s hg18.chr7    27707223 11 + 158545518 agctgaaaaca
        s mm4.chr6     53310104 11 - 151104725 AGCTGAAAATA
        """
        B = """
        a score=2.0 tree='((panTro1.chr6:0.3)baboon:0.2)hg18.chr7;'
        s hg18.chr7    27707221 6 + 158545518 gcagct
        s panTro1.chr6 28869787 5 + 161576975 gcagc-
        s baboon         249182 5 -   4622798 gcagc-
        
        a score=2.0 tree='((panTro1.chr6:0.3)baboon:0.2)hg18.chr7;'
        s hg18.chr7    27707227 7 + 158545518 gaaaaca
        s panTro1.chr6 28869792 7 + 161576975 gaaaaca
        s baboon         249187 7 -   4622798 gaaaaca
        """ 
        C = """
        a score=0.000000 tree='(((panTro1.chr6:0.3)baboon:0.2),mm4.chr6:0.15)hg18.chr7'
        s panTro1.chr6  28869787 12 + 161576975 gcagc-gaaaaca
        s baboon.baboon   249182 12 -   4622798 gcagc-gaaaaca
        s mm4.chr6      53310102 13 - 151104725 ACAGCTGAAAATA
        s hg18.chr7     27707221 13 + 158545518 gcagctgaaaaca
        """
        self.mafJoinTest("hg18", A, B, C)

    def testJoin4(self):
        """Dup join. Contains dups of non-ref sequences and simple split
        """
        A = """
        a score=5.0 tree='((baboon:0.3,baboon:0.1)panTro1.chr6:0.2)hg18.chr7;'
        s baboon         116834 37 +   4622798 AAAGGGAATGTTAACCAAATGAGTTGTCTCTTATGGT
        s baboon         126834 37 +   4622798 AAAGGGAATGTTAACCAAATGAGTTGTCTCTTATGGT
        s panTro1.chr6 28741140 37 + 161576975 AAAGGGAATGTTAACCAAATGAATTGTCTCTTACGGT
        s hg18.chr7    27578828 37 + 158545518 AAAGGGAATGTTAACCAAATGAATTGTCTCTTACGGT
        """
        B = """
        a score=50.0 tree='((mm4.chr6:0.15,mm4.chr7:0.14)rn3.chr4:0.15,rn3.chr4:0.14)hg18.chr7;'
        s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
        s mm4.chr7     50000000 40 + 150000000 AAATGGGAATGTTAAGCAAACGAT--ATTGTCTCTCAGTGTG
        s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
        s rn3.chr4     80000000 35 - 187371129 -AA-GGGGATG-----CCAATGAGTTGTTGTCTCTCAATGTG
        s hg18.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
        """
        C = """
        a score=0.000000 tree='((mm4.chr6:0.15,mm4.chr7:0.14)rn3.chr4:0.15,rn3.chr4:0.14,(baboon:0.3,baboon:0.1)panTro1.chr6:0.2)hg18.chr7'
        s baboon.baboon   126834 37 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGT-
        s baboon.baboon   116834 37 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGT-
        s panTro1.chr6  28741140 37 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGT-
        s mm4.chr6      53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
        s mm4.chr7      50000000 40 + 150000000 AAATGGGAATGTTAAGCAAACGAT--ATTGTCTCTCAGTGTG
        s rn3.chr4      81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
        s rn3.chr4      80000000 35 - 187371129 -AA-GGGGATG-----CCAATGAGTTGTTGTCTCTCAATGTG
        s hg18.chr7     27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
        """
        self.mafJoinTest("hg18", A, B, C)
        
    def testJoin5(self):
        """Dup join contains dups of both ref and non-ref sequences.
        Note the dupped ref human sequences maintain there position in the
        final maf. Once again, the join process should maintain the ordering of sequences in columns.
        This case is interesting as it contains an inserted human base in the dup - I'm allowing
        this, but you might want to consider it illegal and change the test. 
        """
        A = """
        a score=1000.0 tree='((panTro1.chr6:0.3,panTro1.chr6:0.1)baboon:0.2)hg18.chr7;'
        s panTro1.chr6 28862317 6 + 161576975 TAAAGA
        s panTro1.chr6 29862317 6 + 161576975 TAAAGA
        s baboon         241163 3 +   4622798 T---GA
        s hg18.chr7    27699739 3 + 158545518 T---GA
        """
        B = """
        a score=10.0 tree='((hg18.chr7:0.1,hg18.chr7:0.5)mm4.chr6:0.15)rn3.chr4;'
        s hg18.chr7    27699739 3 + 158545518 T---GA
        s hg18.chr7    27000000 5 + 158545518 T-GGGA
        s mm4.chr6     53303881 6 + 151104725 TAAAGA
        s rn3.chr4     81444246 6 + 187371129 taagga
        """
        C = """
        a score=0.000000 tree='((((panTro1.chr6:0.3,panTro1.chr6:0.1)baboon:0.2)hg18.chr7:0.1,hg18.chr7:0.5)mm4.chr6:0.15)rn3.chr4;'
        s panTro1.chr6  29862317 6 + 161576975 T---AAAGA
        s panTro1.chr6  28862317 6 + 161576975 T---AAAGA
        s baboon.baboon   241163 3 +   4622798 T------GA
        s hg18.chr7     27000000 5 + 158545518 T-GG---GA
        s hg18.chr7     27699739 3 + 158545518 T------GA
        s mm4.chr6      53303881 6 + 151104725 TAAA---GA
        s rn3.chr4      81444246 6 + 187371129 taag---ga
        """
        self.mafJoinTest("hg18", A, B, C)
        
    def testJoin6(self):
        """Dup multiple in the reference which joins to sequences..
        """
        A = """
        a score=2.0 tree='((hg18.chr7:0.1,hg18.chr9:0.5)panTro1.chr6:0.15)baboon;'
        s hg18.chr7    27707221 13 + 158545518 gcagctgaaaaca
        s hg18.chr9    27707221 13 + 158545518 gcagctgaaaaca
        s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
        s baboon         249182 13 -   4622798 gcagctgaaaaca
        """
        B = """
        a score=5.0 tree='(mm4.chr6:0.1)hg18.chr7'
        s mm4.chr6     53310102 13 - 151104725 ACAGCTGAAAATA
        s hg18.chr7    27707221 13 + 158545518 gcagctgaaaaca
        
        a score=5.0 tree='(mm4.chr6:0.2)hg18.chr9'
        s mm4.chr6     54310102 13 - 151104725 ACAGCTGAAAATA
        s hg18.chr9    27707221 13 + 158545518 gcagctgaaaaca
        """
        C = """
        a score=0.000000 tree='(((mm4.chr6:0.1)hg18.chr7:0.1,(mm4.chr6:0.2)hg18.chr9:0.5)panTro1.chr6:0.15)baboon;'
        s mm4.chr6      53310102 13 - 151104725 ACAGCTGAAAATA
        s hg18.chr7     27707221 13 + 158545518 gcagctgaaaaca
        s mm4.chr6      54310102 13 - 151104725 ACAGCTGAAAATA
        s hg18.chr9     27707221 13 + 158545518 gcagctgaaaaca
        s panTro1.chr6  28869787 13 + 161576975 gcagctgaaaaca
        s baboon.baboon   249182 13 -   4622798 gcagctgaaaaca
        """
        self.mafJoinTest("hg18", A, B, C)
          
if __name__ == '__main__':
    unittest.main()
