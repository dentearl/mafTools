# mafTools

**mafTools** is a collection of tools that operate on Multiple Alignment Fomat ([maf](http://genome.ucsc.edu/FAQ/FAQformat.html#format5)) files.

## Authors
[Dent Earl](https://github.com/dentearl/), [Benedict Paten](https://github.com/benedictpaten/), [Mark Diekhans](https://github.com/diekhans)

## Dependencies
* sonLib https://github.com/benedictpaten/sonLib/

## Installation
1. Download the package. Consider making it a sibling directory to <code>sonLib/</code>.
2. <code>cd</code> into the directory.
3. Type <code>make</code>.

## Components
* **mafBlockExtractor** A program to extract all alignment blocks that contain a region in a particular sequence. Useful for isolating regions of interest in large maf files.
* **mafBlockFinder** A program to search for a position in a particular sequence. Useful for determining where in maf a particular part of the alignment resides.
* **mafBlockSorter** A program to sort all of the blocks in a MAF based on the (absolute) start position of one of the sequences. Blocks without the sequence are placed at the start of the output in their original order.
* **mafComparator** A program to compare two maf files by sampling. Useful when testing predicted alignments against known true alignments.
* **mafCoveragePickles** A set of programs to assess the pairwise coverage between sequences and to extract the indel distribution of a set of sequences contained in the maf.
* **mafValidator** A program to assess whether or not a given maf file's formatting is valid. 
