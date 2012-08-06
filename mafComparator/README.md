# mafComparator

February 2011 -- August 2012

## Authors

[Dent Earl](https://github.com/dentearl/), [Benedict Paten](https://github.com/benedictpaten/)

## Description
This program takes two [MAF](http://genome.ucsc.edu/FAQ/FAQformat#format5) files and compares them to one another.
Specifically, for each ordered pair of sequences in the first MAF it 
samples a predefined number of sample homology tests (see below), then 
reads the second MAF checking to see which, if any, of the sampled pairs, 
is present. The comparison is then reversed and repeated. Statistics are
then reported in an XML formatted file. MafComparator is suitable for 
running over very large alignments (those with many positions), because 
it does not attempt to hold everything in memory but instead takes a 
sampling approach.

For two sets of pairwise alignments, **A** and **B**, a homology test is 
defined as follows. Pick a pair of aligned positions in **A**, called a 
homology pair -- the **AB** homology test returns _true_ if the pair is present in **B**, 
otherwise it returns _false_. The set of possible homology tests for the 
ordered pair (**A**, **B**) is not necessarily equivalent to the set of 
possible (**B**, **A**) homology tests. We call the proportion of _true_ tests 
(as a percentage of the total of a set of **C** many homology tests), from 
(**A**, **B**) **A~B**.

If **A** is the set of true pairwise alignments and **B** the predicted set of 
alignments then **A~B** (over large enough  **C**), is a proxy to 
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
* <code>-a --logLevel</code> : Set the log level. [off, critical, info, debug] in ascending order
* <code>--maf1</code> : The location of the first MAF file. If comparing true to predicted alignments, this is the truth.
* <code>--maf2</code> : The location of the second MAF file.
* <code>--out</code> : The output XML formatted results file.
* <code>--samples</code> : The ideal number of sample homology tests to perform for the two comparisons (i.e. file1 -> file and file2 -> file1). This number is an ideal because pairs are sampled and thus the actual number may be slightly higher or slightly lower than this value. If this value is equal to or greater than the total number of pairs in a file, then all pairs will be tested. [default 1000000]
* <code>-g --near</code> : The number of bases in either sequence to allow a match to slip by. I.e. <code>--near=n</code> (where _n_ is a non-negative integer) will consider a homology test for a given pair (**S1**:_x_, **S2**:_y_) where **S1** and **S2** are sequences and _x_ and _y_ are positions in the respective sequences, to be a true homology test so long as there is a pair within the other alignment (**S1**:_w_, **S2**:_z_) where EITHER (_w_ is equal to _x_ and _y_ - _n_ <= _z_ <= _y_ + _n_) OR (_x_ - _n_ <= _w_ <= _x_ + _n_ and _y_ is equal to _z_).
* <code>--bedFiles</code> : The location of bed file(s) used to filter the pairwise comparisons. Comma separated list.
* <code>--wigglePairs</code> : The paired names of sequences (comma separated values) to create output that isolates event counts to specific regions of one genome (the first genome in the pair). The asterisk, \*, can be used as wildcard character. i.e. hg19\*,mm9\* will match hg19.chr1 and mm9.chr1 etc etc resulting in all pairs between hg19\* and mm9\*. This feature ignores any intervals described with the <code>--bedFiles</code> option.
* <code>--wiggleBinLength</code> : The length of the bins when the <code>--wigglePairs</code> option is invoked. [default: 100000]
* <code>--numberOfPairs</code> : A pair of comma separated positive integers representing the total number of pairs in maf1 and maf2 (in that order). These numbers are double checked by mafComparator as it runs, a discrpency will cause an error. If these values are known prior to the analysis (either because the analysis has been run before or by use of the mafPairCounter program) this option provides about a 15% speedup. Example: <code>--numberOfPairs 2847390129,228470192212</code>
* <code>--legitSequences</code> : A list of comma separated key value pairs, which themselves are colon (:) separated. Each pair is a sequence name and source length. These values are normally determined by reading all sequences and source lengths from maf1 and then again from maf2 and then finding the intersection of the two sets. The source lengths are verified by mafComparator is it runs and discrepncies will cause errors. If this option is invoked it can result in a speedup of about 15%. Example: <code>--legitSequences apple.chr1:100,apple.chr2:102,pineapple.chr1:2010</code>
* <code>-s --seed</code> : An integer to seed the random number generator. Omitting this causes the seed to be pseudorandom (via <code>time()</code> and <code>getpid()</code>). The seed value is always stored in the output xml.
* <code>-v --version</code> : Print current version number.
* <code>-h --help</code> : Print this help screen.

## Example
Two mafs are included in the example/ directory and can be compared using the command:

    $ mafComparator --maf1 example/a.maf --maf2 example/b.maf --out comparison_a-b.xml

You may note in the output that there are no comparisons for the sequences that are found only in <code>b.maf</code>, i.e. sequences D, E and F. The hash of sequence names used for comparisons is populated using the intersection of the sequence names from the <code>--maf1</code> and <code>--maf2</code> inputs. Sequences that only appear in <code>--maf1</code> or only appear in <code>--maf2</code> input are ignored.
