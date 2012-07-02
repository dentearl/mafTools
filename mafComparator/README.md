# mafComparator

February 2011 -- July 2012

## Authors

[Dent Earl](https://github.com/dentearl/), [Benedict Paten](https://github.com/benedictpaten/)

## Description
This program takes two [MAF](http://genome.ucsc.edu/FAQ/FAQformat#format5) files and for each ordered pair of sequences 
in the MAF it calculates a predefined number of sample homology tests 
(see below), then reports the statistics in an XML formatted file.
It is suitable for running over very large alignments (those with many positions), 
because it does not attempt to hold everything in memory but instead 
takes a sampling approach.

For two sets of pairwise alignments, **A** and **B**, a homology test is 
defined as follows. Pick a pair of aligned positions in **A**, called a 
homology pair -- the **AB** homology test returns _true_ if the pair is present in **B**, 
otherwise it returns _false_. The set of possible homology tests for the 
ordered pair (**A**, **B**) is not necessarily equivalent to the set of 
possible (**B**, **A**) homology tests. We call the proportion of _true_ tests 
(as a percentage of the total of a set of **C** many homology tests), from 
(**A**, **B**) **A~B**.

If **A** is the set of true pairwise alignments and **B** the predicted set of 
alignments then *A~B* (over large enough  **C**), is a proxy to 
[_sensitivity_](http://en.wikipedia.org/wiki/Sensitivity_and_specificity)
of **B** in predicted the set of correctly aligned pairs in **A**. Conversely 
**B~A** (over large enough **C**) is a proxy to the 
[_specificity_](http://en.wikipedia.org/wiki/Sensitivity_and_specificity) of the 
aligned pairs in **B** with respect to the set of correctly aligned pairs 
in **A**.

## Dependencies
* sonLib https://github.com/benedictpaten/sonLib/

## Installation
1. Download the package. Consider making the parent of mafComparator a sibling directory to <code>sonLib</code>.
2. <code>cd</code> into the directory.
3. Type <code>make</code>.

## Use
<code>mafComparator --maf1=FILE1 --maf2=FILE2 --out=OUT.xml [options]</code>

### Options
* <code>mafComparator, version 0.6 July 2012</code>
* <code>-a --logLevel</code> : Set the log level. [off, critical, info, debug] in ascending order.
* <code>--maf1</code> : The location of the first MAF file. If comparing true to predicted alignments, this is the truth.
* <code>--maf2</code> : The location of the second MAF file
* <code>--out</code> : The output XML formatted results file.
* <code>--samples</code> : The number of sample homology tests to perform (total) [default 1000000].
* <code>--bedFiles</code> : The location of bed file(s) used to filter the pairwise comparisons. Comma separated list.
* <code>-g --near</code> : The number of bases in either sequence to allow a match to slip by.
* <code>-s --seed</code> : An integer to seed the random number generator. Omitting this causes the seed to be pseudorandom (via time() and getpid()).
* <code>-v --version</code> : Print current version number
* <code>-h --help</code> : Print this help screen

## Example
Two mafs are included in the example/ directory and can be compared using the command

    $ mafComparator --maf1 example/a.maf --maf2 example/b.maf --out comparison_a-b.xml

You may note in the output that there are no comparisons for the sequences that are found only in <code>b.maf</code>, i.e. sequences D, E and F. The hash of sequence names used for comparisons is populated using the intersection of the sequence names from the <code>--maf1<code> and <code>--maf2<code> inputs. Sequences that only appear in <code>--maf1</code> or only appear in <code>--maf2</code> input are ignored.
