# mafRowOrderer

4 October 2012

## Author

[Dent Earl](https://github.com/dentearl/)

## Description
mafRowOrderer is a program that will look through a maf file block by block and order the maf lines within a block according to the order provided. Species not in the established ordered are excised. Comments are excised. Non sequnece lines ('^s') are excised.

## Installation
1. Download the package.
2. <code>cd</code> into the directory.
3. Type <code>make</code>.

## Use
<code>mafRowOrderer --maf [path to maf] --order [comma separated list of species] </code>

### Options
* <code>-h, --help</code>   show this help message and exit.
* <code>-m, --maf</code>   path to maf file.
* <code>--order</code>   comma separated list of species names
* <code>-v, --verbose</code>   turns on verbose output.

## Example
    $ ./mafRowOrderer --maf example.maf --order hg18,mm9,rn4,banana
    ##maf version=1
    a score=0
    s banana.chr1 0 10 + 1000000 ACGTACGTAC
    ...
    

