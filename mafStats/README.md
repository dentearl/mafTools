# mafStats
5 July 2012

## Author
[Dent Earl](https://github.com/dentearl/)

## Description
A program to read MAF file and report back statistics about the contents.

## Installation
1. Download the package.
2. <code>cd</code> into the directory.
3. Type <code>make</code>.

## Use
<code>mafStats --maf mafFile.maf [options]</code>

### Options
* <code>-h, --help</code>   show this help message and exit.
* <code>-m, --maf</code>     path to maf file.

### Example
    $ mafStats --maf smallDemo.maf
    smallDemo.maf
    ------------------------------
    File size:       66.13 MB
    Lines:             212986
    Header lines:           5
    s lines:           144592
    e lines:                0
    i lines:                0
    q lines:                0
    Blank lines:        68388
    Comment lines:          1
    Sequence chars:     49181016 ( 77.65%)
    Gap chars:          14154166 ( 22.35%)
    Blocks:                      34194
    Ave block area:            1852.23
    Max block area:              37840
    Ave seq field length:       340.14
    Max seq field length:         7568
    Ave seq count in block:       4.23
    Max seq count in block:          5
    10 unique sequences, ordered by # bases present:
                simHuman.chr1:      5311230 ( 10.80%)
                simHuman.chr0:      5225141 ( 10.62%)
                  simDog.chr1:      5048843 ( 10.27%)
                  simDog.chr0:      5023883 ( 10.22%)
                  simCow.chr1:      4989381 ( 10.14%)
                  simCow.chr0:      4979544 ( 10.12%)
                simMouse.chr1:      4671434 (  9.50%)
                simMouse.chr0:      4654920 (  9.46%)
                  simRat.chr1:      4654607 (  9.46%)
                  simRat.chr0:      4622033 (  9.40%)
                        total:     49181016 (100.00%)
    
