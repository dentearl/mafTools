#mafTools

**mafTools** is a collection of tools that operate on .maf files.

##Authors
[Dent Earl](https://github.com/dentearl/), [Benedict Paten](https://github.com/benedictpaten/), Mark Diekhans

##Dependencies
* sonLib https://github.com/benedictpaten/sonLib/
* mafJoin will not compile without the kentSource installed: <code>git clone git://genome-source.cse.ucsc.edu/kent.git</code>

##Installation
1. Download the package. Consider making it a sibling directory to <code>sonLib/</code>.
2. <code>cd</code> into the directory.
3. Type <code>make</code>.
    * If you want to build mafJoin you will need to set kentDir on the command line via <code>make kentDir=/path/to/kent/src</code>
