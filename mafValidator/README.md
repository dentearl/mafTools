#mafValidator

10 October 2011

##Authors

[Dent Earl](https://github.com/dentearl/)

##Description
mafValidator is a script to validate the formatting and basic data contained in a maf file.

##Dependencies
* Python 2.6 &le; version &lt; 3.0

##Installation
1. Download the package.
2. <code>cd</code> into the directory.
3. Type <code>make</code>.

##Use
<code>mafValidator.py --maf=FILE [options]</code>

###Options
* <code>mafValidator.py</code>
* <code>--maf</code> : path to the maf file to test
* <code>--testChromNames</code> : Expects that the source field will be formatted with .chrN e.g. "hg19.chr1"
* <code>--help</code> : Print this help screen

##Test
<code>make test</code>
