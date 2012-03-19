#mafBlockSorter

19 March 2012

##Author

[Dent Earl](https://github.com/dentearl/)

##Description
mafBlockSorter is a program that will sort the blocks of a maf in ascending order of the sequence start field of the specified sequence name. Blocks that do not contain the specified sequence will be output at the start of the maf in the order they appear in the input, followed by the sorted blocks.

##Dependencies
* 

##Installation
1. Download the package.
2. <code>cd</code> into the directory.
3. Type <code>make</code>.

##Use
<code>mafBlockSorter --seq [sequence name (and possibly chr)] [options] < myFile.maf</code>

###Options
* <code>-h, --help</code>   show this help message and exit.
* <code>-s, --seq</code>   sequence _name.chr_ e.g. `hg18.chr2'.
* <code>-v, --verbose</code>   turns on verbose output.

##Example
    $ ./mafBlockSorter --seq hg19.chr20 < example.maf 
    ##maf version=1 
    
    #a score=0 pctid=99.2
    #s hg19.chr20 0 795 + 73767698 GAT...
    ...

