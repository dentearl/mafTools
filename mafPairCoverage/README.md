# mafPairCoverage

7 May 2013

## Author

[Dent Earl](https://github.com/dentearl/)

## Description
mafPairCoverage is a program that will look through a maf file block by block and check for a particular pair of sequences (allowing input sequence to end in wildcard *) and count the number of aligned positions where the two sequences have residues aligned. Coverage of genome A onto genome B is then symmetrically calculated as the number of aligned positions divided by the total size of genome B.

__BE AWARE!__ The input maf should be transitively closed (if you are unsure you can use the tool mafTransitiveClosure to transitively close the alignment) to insure that the coverage numbers are accurate.

## Installation
1. Download the package.
2. <code>cd</code> into the directory.
3. Type <code>make</code>.

## Use
<code>mafPairCoverage --seq1 [sequence name] --seq2 [sequence name] --maf myFile.maf [options]</code>

### Options
* <code>-h, --help</code>   show this help message and exit.
* <code>--seq1</code>   sequence _name.chr_ e.g. `hg19*'. May end in * to indicate wildcard.
* <code>--seq2</code>   sequence _name.chr_ e.g. `mm9.chr2'. May end in * to indicate wildcard.
* <code>--maf</code>    input maf file.
* <code>-v, --verbose</code>   turns on verbose output.

## Example
    $ ./mafPairCoverage --seq1 hg19* --seq2 mm9* --maf example.maf 
    ...

