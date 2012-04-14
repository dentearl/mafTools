#mafComparator

February -- November 2011

##Authors

[Benedict Paten](https://github.com/benedictpaten/), [Dent Earl](https://github.com/dentearl/)

##Description
The script takes two MAF files and for each ordered pair of sequences 
in the MAFS calculates a predefined number of sample homology tests 
(see below), then reports the statistics in an XML formatted file.
It is suitable for running over very large numbers of alignments, 
because it does not attempt to hold everything in memory, and instead 
takes a sampling approach.

For two sets of pairwise alignments, A and B, a homology test is 
defined as follows. Pick a pair of aligned positions in A, called a 
homology pair - the AB homology test returns true if the pair is in B, 
otherwise it returns false. The set of possible homology tests for the 
ordered pair (A, B) is not necessarily equivalent to the set of 
possible (B, A) homology tests. We call the proportion of true tests 
as a percentage of the total of a set of homology tests C from 
(A, B)  A~B.

If A is the set of true pairwise alignments and B the predicted set of 
alignments then A~B (over large enough  C), is a proxy to sensitivity 
of B in predicted the set of correctly aligned pairs in A. Conversely 
B~A (over large enough C) is a proxy to the specificity of the 
aligned pairs in B with respect to the set of correctly aligned pairs 
in A.

##Dependencies
* sonLib https://github.com/benedictpaten/sonLib/

##Installation
1. Download the package. Consider making the parent of mafComparator a sibling directory to <code>sonLib</code>.
2. <code>cd</code> into the directory.
3. Type <code>make</code>.

##Use
<code>mafComparator --mafFile1=FILE1 --mafFile2=FILE2 --outputFile=OUT.xml [options]</code>

###Options
* <code>mafComparator, version 0.3</code>
* <code>-a --logLevel</code> : Set the log level. [off, critical, info, debug] in ascending order.
* <code>-b --mafFile1</code> : The location of the first MAF file. If comparing true to predicted alignments, this is the truth.
* <code>-c --mafFile2</code> : The location of the second MAF file
* <code>-d --outputFile</code> : The output XML formatted results file.
* <code>-e --sampleNumber</code> : The number of sample homology tests to perform (total) [default 1000000].
* <code>-p --printFailures</code> : Print tab-delimited details about failed tests to stderr.
* <code>-f --bedFiles</code> : The location of bed file(s) used to filter the pairwise comparisons. Comma separated list.
* <code>-g --near</code> : The number of bases in either sequence to allow a match to slip by.
* <code>-s --seed</code> : An integer to seed the random number generator. Omitting this causes the seed to be pseudorandom (via time() and getpid()).
* <code>-v --version</code> : Print current version number
* <code>-h --help</code> : Print this help screen

##Example
Two mafs are included in the example/ directory and can be compared using the command

<code>$ mafComparator --mafFile1 example/a.maf --mafFile2 example/b.maf --outputFile comparison_a-b.xml</code>

You may note in the output that there are no comparisons for the sequences that are found only in b.maf, i.e. sequences D, E and F. The hash of sequence names used for comparisons is populated using the intersection of the sequence names from the --mafFile1 and --mafFile2 inputs. Sequences that only appear in --mafFile1 or only appear in --mafFile2 input are ignored.
