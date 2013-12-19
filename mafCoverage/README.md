# mafCoverage

December 2013

## Author

[Benedict Paten](https://github.com/benedictpaten/)

## Description
mafCoverage is a program that will look through a maf file block by block and check for the coverage of all other sequences onto one user-specified sequence.

The input need not be transitively closed as <code>mafCoverage</code> builds a bit array for the user-specied sequence and stores only presence-absense data. Duplications are only counted once.

## Installation
1. Download the package.
2. <code>cd</code> into the directory.
3. Type <code>make</code>.

## Use
<code>mafCoverage ... file this in</code>

### Options
* <code>-h, --help</code>   show this help message and exit.
File this in

## Example
    $ ./mafCoverage ... file this in
    ...

