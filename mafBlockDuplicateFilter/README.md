# mafBlockDuplicateFilter

16 April 2012

## Author
[Dent Earl](http://github.com/dentearl/)

## Description
mafBlockDuplicateFilter is a program to filter out duplications from a Multiple Alignment Format (maf) file. This program assumes the sequence name field is formatted as in "speciesName.chromosomeName" using the first period charater, ".", as the delimiter between the species name and the chromosome name. For every block present in the alignment, mBDF looks for any duplicated species within the block. Instead of stripping out all copies of the duplication, the sequence with the highest similarity to the consensus of the block is left, all others are removed. Sequence similarity is computed as a bit score in comparison to the IUPAC-enabled consensus. Ties are resolved by picking the sequence that appears earliest in the file.

## Installation
1. Download the package.
2. <code>cd</code> into the directory.
3. Type <code>make</code>.

## Use
<code>mafBlockDuplicateFilter < mafWithDuplicates.maf > pruned.maf </code>

### Options
* <code>-h, --help</code>   show this help message and exit.

## Example
    $ ./mafBlockDuplicateFilter < mafWithDuplicates.maf > mafPruned.maf

