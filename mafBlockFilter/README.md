# mafBlockFilter

28 May 2012

## Author

[Dent Earl](https://github.com/dentearl/)

## Description
mafBlockFilter is a program that will look through a maf file block by block and excise out sequence lines that match criteria established by the user on the command line. For example one can filter out all sequence lines that start with 'hg18' using <code>--exclude</code> or filter for sequence lines starting with only 'hg19', 'mm9' and 'rn4' using <code>--include</code>.

## Installation
1. Download the package.
2. <code>cd</code> into the directory.
3. Type <code>make</code>.

## Use
<code>mafBlockFilter --maf [path to maf] [options] </code>

### Options
* <code>-h, --help</code>   show this help message and exit.
* <code>-m, --maf</code>   path to maf file.
* <code>-i, --includeSeq</code>   comma separated list of sequence names to include
* <code>-e, --excludeSeq</code>   comma separated list of sequence names to exclude
* <code>-g, --noDegreeGT</code>       filter out all blocks with degree greater than this value.
* <code>-l, --noDegreeLT</code>       filter out all blocks with degree less than this value.
* <code>-v, --verbose</code>   turns on verbose output.

## Example
    $ ./mafBlockFilter --maf example.maf --include hg18,mm9,rn4,banana
    ##maf version=1
    a score=0
    s banana.chr1 0 10 + 1000000 ACGTACGTAC
    ...
    

