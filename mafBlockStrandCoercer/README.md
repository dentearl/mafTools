# mafBlockStrandCoercer

3 October 2012

## Author
[Dent Earl](http://github.com/dentearl/)

## Description
mafBlockStrandCoercer is a program to coerce a particular strandedness out for a block based the strandedness of a target sequence. When a block contains the target sequence but in the flipped orientation then the block is flipped, i.e. all start coordinates are transformed, and all sequence fields are reverse-complemented. If the block contains the target sequence multiple times and with conflicing strands (i.e. both + and - strands are observed), then nothing is done.

## Installation
1. Download the package.
2. <code>cd</code> into the directory.
3. Type <code>make</code>.

## Use
<code>mafBlockStrandCoercer --maf alignment.maf --seq hg18 --strand + > positive.maf </code>

### Options
* <code>-h, --help</code>   show this help message and exit.
* <code>--maf</code>   input alignment maf file.
* <code>--seq</code>   sequence to base block strandedness upon. (string comparison only done for length of input, i.e. --seq=hg18 will match hg18.chr1, hg18.chr2, etc etc)
* <code>--strand</code>   strand to enforce, when possible. may be + or -, defaults to +.

## Example
    $ mafBlockStrandCoercer --maf alignment.maf --seq hg18 --strand + > positive.maf 

