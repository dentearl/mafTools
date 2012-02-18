#mafBlockExtractor

14 Feb 2012

##Author

[Dent Earl](https://github.com/dentearl/)

##Description
mafBlockExtractor is a program that will look through a maf file for a particular sequence name and region. If a match is found then the block containing the querry will be printed to standard out.

##Dependencies
* 

##Installation
1. Download the package.
2. <code>cd</code> into the directory.
3. Type <code>make</code>.

##Use
<code>mafBlockFinder --seq [sequence name (and possibly chr)] --pos [position to search for] [options] < myFile.maf</code>

###Options
* <code>-h, --help</code>   show this help message and exit.
* <code>-s, --seq</code>   sequence _name.chr_ e.g. `hg18.chr2'.
* <code>--start</code>   start of the region, inclusive. Must be a positive number.
* <code>--stop</code>   end of the region, inclusive. Must be a positive number.
* <code>-v, --verbose</code>   turns on verbose output.

##Example
    $ ./mafBlockExractor --seq hg19.chr20 --start 500 --stop 1000 < example.maf 
    ##maf version=1 
    
    #a score=0 pctid=99.2
    #s hg19.chr20 0 795 + 73767698 GAT...
    ...

