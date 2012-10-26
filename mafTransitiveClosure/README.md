# mafTransitiveClosure
24 May 2012

## Author
[Dent Earl](https://github.com/dentearl/)

## Description
A program to perform the transitive closure on an alignment. That is it checks every column of the alignment and looks for situations where a position A is aligned to B in one part of a file and B is aligned to C in another part of the file. The transitive closure of this relationship would be a single column with A, B and C all present. Useful for when you have pairwise alignments and you wish to turn them into something more resembling a multiple alignment.

## Dependencies
* sonLib https://github.com/benedictpaten/sonLib/
* pinchesAndCacti https://github.com/benedictpaten/pinchesAndCacti

## Installation
1. Download the package.
2. <code>cd</code> into the directory.
3. Type <code>make</code>.

## Use
<code>mafTransitiveClosure --maf mafFile.maf > transitivelyClosed.maf</code>

### Options
* <code>-h, --help</code>   show this help message and exit.
* <code>-m, --maf</code>     path to maf file.
* <code>-v, --verbose</code>   turns on verbose output.
