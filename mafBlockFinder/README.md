# mafBlockFinder

10 Feb 2012

## Author

[Dent Earl](https://github.com/dentearl/)

## Description
mafBlockFinder is a program that will look through a maf file for a particular sequence name and location. If a match is found the line number and first few fields are returned. If no match is found nothing is returned.

## Installation
1. Download the package.
2. <code>cd</code> into the directory.
3. Type <code>make</code>.

## Use
<code>mafBlockFinder --maf [path to maf] --seq [sequence name (and possibly chr)] --pos [position to search for, zero based coordinates] [options] </code>

### Options
* <code>-h, --help</code>   show this help message and exit.
* <code>-m, --maf</code>   path to maf file.
* <code>-s, --seq</code>   sequence _name.chr_ e.g. `hg18.chr2'.
* <code>-p, --pos</code>   position along the chromosome you are searching for. Must be a non negative number.
* <code>-v, --verbose</code>   turns on verbose output.

## Example
    $ ./mafBlockFinder --maf example.maf --seq apple.chr20 --pos 500
    4: s apple.chr20 0 795 + 73767698 ...AATTG ->G<- ACCCG...
    
We see from this example that position 500 of apple.chr20 is located at line 4 of example.maf and that the base at this position is G flanked by AATTG on the left and ACCCG on the right.
