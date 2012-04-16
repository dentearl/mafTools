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
<code>mafBlockFinder --seq [sequence name (and possibly chr)] --pos [position to search for] [options] < myFile.maf</code>

### Options
* <code>-h, --help</code>   show this help message and exit.
* <code>-s, --seq</code>   sequence _name.chr_ e.g. `hg18.chr2'.
* <code>-p, --pos</code>   position along the chromosome you are searching for. Must be a positive number.
* <code>-v, --verbose</code>   turns on verbose output.

## Example
    $ ./mafBlockFinder --seq hg19.chr20 --pos 500 < example.maf 
    4: s hg19.chr20 0 795 + 73767698 ...

