#mafComparator

16 February 2011

##Authors

[Benedict Paten](https://github.com/benedictpaten/), [Dent Earl](https://github.com/dentearl/)

##Dependencies
* sonLib https://github.com/benedictpaten/sonLib/

##Installation
1. Download the package. Consider making the parent of mafComparator a sibling directory to <code>sonLib</code>.
2. <code>cd</code> into the directory.
3. Type <code>make</code>.

##Use
<code>mafComparator --mafFile1=FILE1 --mafFile2=FILE2 --outputFile=OUT.xml [options]</code>

###Options
* <code>mafComparator, version 0.1</code>
* <code>-a --logLevel</code> : Set the log level. [off, critical, info, debug] in ascending order.
* <code>-b --mafFile1</code> : The location of the first MAF file (used to create sequence name hash.)
* <code>-c --mafFile2</code> : The location of the second MAF file
* <code>-d --outputFile</code> : The output XML formatted results file.
* <code>-e --sampleNumber</code> : The number of sample homology tests to perform (total) [default 1000000].
* <code>-p --printFailures</code> : Print tab-delimited details about failed tests to stderr.
* <code>-v --version</code> : Print current version number
* <code>-h --help</code> : Print this help screen
* <code>-f --bedFiles</code> : The location of bed file used to filter the pairwise comparisons.
* <code>-g --near</code> : The number of bases in either sequence to allow a match to slip by.

##Example
Two mafs are included in the example/ directory and can be compared using the command

<code>$ mafComparator --mafFile1 example/a.maf --mafFile2 example/b.maf --outputFile comparison_a-b.xml</code>

You may note in the output that there are no comparisons for the sequences that are found only in b.maf, i.e. sequences D, E and F. The hash of sequence names used for comparisons is populated using the --mafFile1 input. Sequences that only appear in the --mafFile2 input are ignored.
