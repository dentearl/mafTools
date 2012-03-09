import mafValidator as mafval
import os
import unittest
import shutil

goodMafs = ['''##maf version=1 scoring=tba.v8
# tba.v8 (((human chimp) baboon) (mouse rat))

a score=23262.0
s hg16.chr7    27707221 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon.chr0    249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA

a score=0
s hg16.chr7    27707221 13 + 158545518 gcagctgaaaaca 
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
i panTro1.chr6 N 0 C 0
s baboon.chr0    249182 13 +   4622798 gcagctgaaaaca
i baboon.chr0   I 234 n 19

a score=0
s hg16.chr7    27707221 13 + 158545518 gcagctgaaaaca
e mm4.chr6     53310102 13 + 151104725 I

a score=0
s hg18.chr1                  32741 26 + 247249719 TTTTTGAAAAACAAACAACAAGTTGG
s panTro2.chrUn            9697231 26 +  58616431 TTTTTGAAAAACAAACAACAAGTTGG
q panTro2.chrUn                                   99999999999999999999999999
s dasNov1.scaffold_179265     1474  7 +      4584 TT----------AAGCA---------
q dasNov1.scaffold_179265                         99----------32239--------- 

''',
            '''##maf version=1 scoring=tba.v8 
# tba.v8 (((human chimp) baboon) (mouse rat)) 
                   
a score=23262.0     
s hg18.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon.chr0    116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
                   
a score=5062.0                    
s hg18.chr7    27699739 6 + 158545518 TAAAGA
s panTro1.chr6 28862317 6 + 161576975 TAAAGA
s baboon.chr0    241163 6 +   4622798 TAAAGA 
s mm4.chr6     53303881 6 + 151104725 TAAAGA
s rn3.chr4     81444246 6 + 187371129 taagga

a score=6636.0
s hg18.chr7    27707221 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon.chr0    249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA

''']

def testFile(s):
   makeTempDir()
   mafFile = os.path.join(os.curdir, 'tempTestDir', 'test.maf')
   f = open(mafFile, 'w')
   f.write(s)
   f.close()
   return mafFile

def makeTempDir():
   if not os.path.exists(os.path.join(os.curdir, 'tempTestDir')):
      os.mkdir(os.path.join(os.curdir, 'tempTestDir'))

def removeTempDir():
   if os.path.exists(os.path.join(os.curdir, 'tempTestDir')):
      shutil.rmtree(os.path.join(os.curdir, 'tempTestDir'))

class HeaderCheck(unittest.TestCase):
   badHeaders = ['''#maf version=1 scoring=tba.v8 
# tba.v8 (((human chimp) baboon) (mouse rat)) 
                   
a score=23262.0     
s hg18.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon.chr0    116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG

''',
                 '''maf version=1 scoring=tba.v8 
                   
a score=23262.0     
s hg18.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon.chr0    116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG

''',
                 '''##maf version = 1 scoring = tba.v8 

a score=23262.0     
s hg18.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon.chr0    116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG

''',
                 '''##maf scoring=tba.v8 

a score=23262.0     
s hg18.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon.chr0    116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG

''',
                 ]
         
   def testBadHeaders(self):
      """ mafValidator should fail on bad headers
      """
      for b in self.badHeaders:
         makeTempDir()
         mafFile = testFile(b)
         self.assertRaises(mafval.HeaderError, mafval.validateMaf, mafFile)
         removeTempDir()
      

class FooterCheck(unittest.TestCase):
   badFooters = ['''##maf version=1 scoring=tba.v8 

a score=23262.0     
s hg18.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon.chr0    116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
''',
                 '''##maf version=1 scoring=tba.v8 

a score=23262.0     
s hg18.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon.chr0    116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG''',
      ]
   def testBadFooters(self):
      """ mafValidator should fail on bad footers
      """
      for b in self.badFooters:
         makeTempDir()
         mafFile = testFile(b)
         self.assertRaises(mafval.FooterError, mafval.validateMaf, mafFile)
         removeTempDir()

class FieldNumberCheck(unittest.TestCase):
   badFieldNumbers = ['''##maf version=1 scoring=tba.v8 

a score=23262.0     
s hg18.7       27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon.chr0    116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344    + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG

''',
                 '''##maf version=1 scoring=tba.v8 

a score=23262.0     
s hg18         27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG      ATTT
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon.chr0    116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG

''',
      ]
   def testFieldNumbers(self):
      """mafValidator should fail when the number of fields in a seq line is not 7
      """
      for b in self.badFieldNumbers:
         mafFile = testFile(b)
         self.assertRaises(mafval.FieldNumberError, mafval.validateMaf, mafFile)
         removeTempDir()

class StrandCheck(unittest.TestCase):
   badStrands = ['''##maf version=1 scoring=tba.v8 

a score=23262.0     
s hg18.chr1    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon.chr0    116834 38 -1  4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG

''',
                 '''##maf version=1 scoring=tba.v8 

a score=23262.0     
s hg18.chr1    27578828 38 1 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon.chr0    116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG

''',
      ]
   def testStrands(self):
      """mafValidator should fail when strand field is malformed
      """
      for b in self.badStrands:
         mafFile = testFile(b)
         self.assertRaises(mafval.StrandCharacterError, mafval.validateMaf, mafFile)
         removeTempDir()

class NamesCheck(unittest.TestCase):
   badNames = ['''##maf version=1 scoring=tba.v8 

a score=23262.0     
s hg18.7       27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon.chr0    116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG

''',
                 '''##maf version=1 scoring=tba.v8 

a score=23262.0     
s hg18         27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon.chr0    116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG

''',
      ]
   def testSrcNames(self):
      """mafValidator should fail when name (src) field is malformed
      """
      for b in self.badNames:
         mafFile = testFile(b)
         self.assertRaises(mafval.SpeciesFieldError, mafval.validateMaf, mafFile, testChromNames = True)
         removeTempDir()

class SourceLengthChecks(unittest.TestCase):
   badSources = ['''##maf version=1 scoring=tba.v8
# tba.v8 (((human chimp) baboon) (mouse rat))

a score=23262.0
s hg16.chr7    27707221 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon.chr0    249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA

a score=23262.0
s hg16.chr7           0 13 +       100 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon.chr0    249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA

''',
                 ]

   def testSourceLengths(self):
      """mafValidator should fail when source length changes
      """
      for i, b in enumerate(self.badSources):
         mafFile = testFile(b)
         if b == 0:
            self.assertRaises(mafval.SourceLengthError, mafval.validateMaf, mafFile, testChromNames = True)
         self.assertRaises(mafval.SourceLengthError, mafval.validateMaf, mafFile)
         removeTempDir()

class AlignmentLengthChecks(unittest.TestCase):
   badLengths = ['''##maf version=1 scoring=tba.v8 

a score=23262.0     
s hg18.7       27578828 39 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon.chr0    116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG

''',
                 '''##maf version=1 scoring=tba.v8 

a score=23262.0     
s hg18         27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon.chr0    116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTGG

''',
      ]
   def testAlignmentLengths(self):
      """mafValidator should fail when seq len is not equal to non dash characters in alignment field
      """
      for b in self.badLengths:
         mafFile = testFile(b)
         self.assertRaises(mafval.AlignmentLengthError, mafval.validateMaf, mafFile)
         removeTempDir()

class AlignmentFieldLengthChecks(unittest.TestCase):
   badLengths = ['''##maf version=1 scoring=tba.v8 

a score=23262.0     
s hg18.7       27578828 42 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTGacgt
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon.chr0    116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG

''',
                 '''##maf version=1 scoring=tba.v8 

a score=23262.0     
s hg18         27578828 36 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon.chr0    116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 42 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTGGG

''',
      ]
   def testAlignmentLengths(self):
      """mafValidator should fail when alignment fields in a block are not all of equal length
      """
      for b in self.badLengths:
         mafFile = testFile(b)
         self.assertRaises(mafval.AlignmentLengthError, mafval.validateMaf, mafFile)
         removeTempDir()

class AlignmentBlockLinesChecks(unittest.TestCase):
   badAlignmentBlockStarts = ['''##maf version=1 scoring=tba.v8 

a score=23262.0
s hg18.7       27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon.chr0    116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG

s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG

''',
                 '''##maf version=1 scoring=tba.v8 

a score=23262.0
s hg18         27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon.chr0    116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG

a score=23262.0
s hg18         27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG

s hg18         27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG

''',
      ]
   badAlignmentBlockKeyValuePairs = ['''##maf version=1 scoring=tba.v8 

a score=23262.0
s hg18.7       27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon.chr0    116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG

a score=0 banana
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 39 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG

''',
                 '''##maf version=1 scoring=tba.v8 

a score=23262.0
s hg18         27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon.chr0    116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG

a score=23262.0
s hg18         27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG

a score=0 banana=apple cheeseburger
s hg18         27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG

''',
      ]
   def testAlignmentBlockLineExistence(self):
      """mafValidator should fail when a sequence block starts without an '^a' line
      """
      for b in self.badAlignmentBlockStarts:
         mafFile = testFile(b)
         self.assertRaises(mafval.MissingAlignmentBlockLineError, mafval.validateMaf, mafFile)
         removeTempDir()
   def testAlignmentBlockLineKeyValuePairs(self):
      """mafValidator should fail when an alignment block has mal-formed key-value pairs
      """
      for b in self.badAlignmentBlockKeyValuePairs:
         mafFile = testFile(b)
         self.assertRaises(mafval.AlignmentBlockLineKeyValuePairError, mafval.validateMaf, mafFile)
         removeTempDir()
   

class StartFieldChecks(unittest.TestCase):
   badFields = ['''##maf version=1 scoring=tba.v8 
# start field cannot be negative

a score=23262.0     
s hg18         27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon.chr0    116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
s dae.chr0         -100 30 +       100 AAA-GGGAATGTTAACCAAATGA-----------TTACGGTG

''',]
   def testFields(self):
      """mafValidator should fail when start fields contain negative values
      """
      for b in self.badFields:
         mafFile = testFile(b)
         self.assertRaises(mafval.StartFieldError, mafval.validateMaf, mafFile)
         removeTempDir()

class SourceSizeFieldChecks(unittest.TestCase):
   badFields = ['''##maf version=1 scoring=tba.v8 
# src length cannot be negative

a score=23262.0     
s hg18         27578828  38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140  38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon.chr0    116834  38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344  38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243  40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
s dae.chr0          100  30 +      -100 AAA-GGGAATGTTAACCAAATGA-----------TTACGGTG

''',]
   def testFields(self):
      """mafValidator should fail when source size fields contain negative values
      """
      for b in self.badFields:
         mafFile = testFile(b)
         self.assertRaises(mafval.SourceSizeFieldError, mafval.validateMaf, mafFile)
         removeTempDir()

class RangeChecks(unittest.TestCase):
   badRanges = ['''##maf version=1 scoring=tba.v8 

a score=23262.0     
s hg18.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon.chr0    116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
s dae.chr0            0 11 +        10 A----------CTAAGCCAA---------------------G

''',
                 '''##maf version=1 scoring=tba.v8 

a score=23262.0     
s hg18         27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon.chr0    116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
s dae.chr0            0 11 -       10  A----------CTAAGCCAA---------------------G

''',
                '''##maf version=1 scoring=tba.v8 

a score=23262.0     
s hg18         27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon.chr0    116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
s dae.chr0            5  6 -        10 A----------CTAA--------------------------G

''',
                ]
   def testSequenceRanges(self):
      """mafValidator should fail when sequences ranges go outside of source length
      """
      for b in self.badRanges:
         mafFile = testFile(b)
         self.assertRaises(mafval.OutOfRangeError, mafval.validateMaf, mafFile)
         removeTempDir()

class LinesStartingWithIChecks(unittest.TestCase):
   badBlocks = ['''##maf version=1 scoring=tba.v8 
# wrong number of fields on i line

a score=23262.0     
s hg18         27578828  38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140  38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
i panTro1.chr6 N 0 C 0 extra
s baboon.chr0    116834  38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344  38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243  40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
s dae.chr0            0  30 +       100 AAA-GGGAATGTTAACCAAATGA-----------TTACGGTG

''',
                '''##maf version=1 scoring=tba.v8 
# non numerics for fields 4 and 6 on the i line

a score=23262.0     
s hg18         27578828  38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140  38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
i panTro1.chr6 N 0 C steve
s baboon.chr0    116834  38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344  38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243  40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
s dae.chr0            0  30 +       100 AAA-GGGAATGTTAACCAAATGA-----------TTACGGTG

''',
                '''##maf version=1 scoring=tba.v8 
# non numerics for fields 4 and 6 on the i line

a score=23262.0     
s hg18         27578828  38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140  38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
i panTro1.chr6 N brian C steve
s baboon.chr0    116834  38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344  38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243  40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
s dae.chr0            0  30 +       100 AAA-GGGAATGTTAACCAAATGA-----------TTACGGTG

''',
                '''##maf version=1 scoring=tba.v8 
# non numerics for fields 4 and 6 on the i line

a score=23262.0     
s hg18         27578828  38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140  38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
i panTro1.chr6 N william C 0
s baboon.chr0    116834  38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344  38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243  40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
s dae.chr0            0  30 +       100 AAA-GGGAATGTTAACCAAATGA-----------TTACGGTG

''',
                '''##maf version=1 scoring=tba.v8 
# negative numbers in the i line numerical fields

a score=23262.0     
s hg18         27578828  38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140  38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
i panTro1.chr6 N -1 C -5
s baboon.chr0    116834  38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344  38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243  40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
s dae.chr0            0  30 +       100 AAA-GGGAATGTTAACCAAATGA-----------TTACGGTG

''',
                '''##maf version=1 scoring=tba.v8 
# negative numbers in the i line numerical fields

a score=23262.0     
s hg18         27578828  38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140  38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
i panTro1.chr6 N -1 C -5
s baboon.chr0    116834  38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344  38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243  40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
s dae.chr0            0  30 +       100 AAA-GGGAATGTTAACCAAATGA-----------TTACGGTG

''',
                '''##maf version=1 scoring=tba.v8 
# negative numbers in the i line numerical fields

a score=23262.0     
s hg18         27578828  38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140  38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
i panTro1.chr6 N 0 C -5
s baboon.chr0    116834  38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344  38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243  40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
s dae.chr0            0  30 +       100 AAA-GGGAATGTTAACCAAATGA-----------TTACGGTG

''',
                '''##maf version=1 scoring=tba.v8 
# negative numbers in the i line numerical fields

a score=23262.0     
s hg18         27578828  38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140  38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
i panTro1.chr6 N -1 C 0
s baboon.chr0    116834  38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344  38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243  40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
s dae.chr0            0  30 +       100 AAA-GGGAATGTTAACCAAATGA-----------TTACGGTG

''',
                '''##maf version=1 scoring=tba.v8 
# floats in the i line numerical fields

a score=23262.0     
s hg18         27578828  38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140  38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
i panTro1.chr6 N 0.7 C 0.35
s baboon.chr0    116834  38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344  38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243  40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
s dae.chr0            0  30 +       100 AAA-GGGAATGTTAACCAAATGA-----------TTACGGTG

''',
                '''##maf version=1 scoring=tba.v8 
# floats in the i line numerical fields

a score=23262.0     
s hg18         27578828  38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140  38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
i panTro1.chr6 I 255.7 I 9.35
s baboon.chr0    116834  38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344  38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243  40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
s dae.chr0            0  30 +       100 AAA-GGGAATGTTAACCAAATGA-----------TTACGGTG

''',
                '''##maf version=1 scoring=tba.v8 
# floats in the i line numerical fields

a score=23262.0     
s hg18         27578828  38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140  38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
i panTro1.chr6 I 255 I 9.35
s baboon.chr0    116834  38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344  38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243  40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
s dae.chr0            0  30 +      100 AAA-GGGAATGTTAACCAAATGA-----------TTACGGTG

''',
                '''##maf version=1 scoring=tba.v8 
# floats in the i line numerical fields

a score=23262.0     
s hg18         27578828  38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140  38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
i panTro1.chr6 I 255.000001 I 9
s baboon.chr0    116834  38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344  38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243  40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
s dae.chr0            0  30 +       100 AAA-GGGAATGTTAACCAAATGA-----------TTACGGTG

''',
                '''##maf version=1 scoring=tba.v8 
# invalid Status codes for i line

a score=23262.0     
s hg18         27578828  38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140  38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
i panTro1.chr6 & 255 C 0
s baboon.chr0    116834  38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344  38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243  40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
s dae.chr0            0  30 +      100 AAA-GGGAATGTTAACCAAATGA-----------TTACGGTG

''',
                '''##maf version=1 scoring=tba.v8 
# invalid Status codes for i line

a score=23262.0     
s hg18         27578828  38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140  38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
i panTro1.chr6 I 255 z 0
s baboon.chr0    116834  38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344  38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243  40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
s dae.chr0            0  30 +      100 AAA-GGGAATGTTAACCAAATGA-----------TTACGGTG

''',
                '''##maf version=1 scoring=tba.v8 
# invalid Status codes for i line

a score=23262.0     
s hg18         27578828  38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140  38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
i panTro1.chr6 a 255 M 0
s baboon.chr0    116834  38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344  38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243  40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
s dae.chr0            0  30 +       100 AAA-GGGAATGTTAACCAAATGA-----------TTACGGTG

''',
                '''##maf version=1 scoring=tba.v8 
# I status must yield value

a score=23262.0     
s hg18         27578828  38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140  38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
i panTro1.chr6 I 0 C 0
s baboon.chr0    116834  38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344  38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243  40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
s dae.chr0            0  30 +       100 AAA-GGGAATGTTAACCAAATGA-----------TTACGGTG

''',
                '''##maf version=1 scoring=tba.v8 
# i line has wrong src field name

a score=23262.0     
s hg18         27578828  38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140  38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
i banana I 1 C 0
s baboon.chr0    116834  38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344  38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243  40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
s dae.chr0            0  30 +       100 AAA-GGGAATGTTAACCAAATGA-----------TTACGGTG

''',
                '''##maf version=1 scoring=tba.v8 
# i line must follow s line

a score=23262.0
i banana I 1 C 0
s hg18         27578828  38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140  38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon.chr0    116834  38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344  38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243  40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
s dae.chr0            0  30 +       100 AAA-GGGAATGTTAACCAAATGA-----------TTACGGTG

''',
                ]
   def testFields(self):
      """mafValidator should fail when "i" lines are malformed
      """
      for b in self.badBlocks:
         mafFile = testFile(b)
         self.assertRaises(mafval.ILineFormatError, mafval.validateMaf, mafFile)
         removeTempDir()

class GoodKnownMafs(unittest.TestCase):
   def testGoodMafs(self):
      """ mafValidator should accept known good mafs
      """
      for g in goodMafs:
         makeTempDir()
         mafFile = testFile(g)
         self.assertTrue(mafval.validateMaf(mafFile))
         removeTempDir()

if __name__ == '__main__':
   unittest.main()
