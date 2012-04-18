# mafBlockDuplicateFilter

16 April 2012

## Author
[Dent Earl](http://github.com/dentearl/)

## Description
mafBlockDuplicateFilter is a program to filter out duplications from a Multiple Alignment Format (maf) file. This program assumes the sequence name field is formatted as in "speciesName.chromosomeName" using the first period charater, ".", as the delimiter between the species name and the chromosome name. For every block present in the alignment, mBDF looks for any duplicated species. Instead of stripping out all copies of the duplication, the sequence with the highest similarity to the consensus of the block is left but all others are removed.

## Installation
1. Download the package.
2. <code>cd</code> into the directory.
3. Type <code>make</code>.

## Use
<code></code>

### Options
* <code>-h, --help</code>   show this help message and exit.
* <code></code>
* <code></code>

## Example
    $ ./mafBlockDuplicateFilter < mafWithDuplicates.maf > mafPruned.maf

