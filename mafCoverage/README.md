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
```shell
Usage: mafCoverage [maf file]

Reports the pairwise (n-)coverage between a specified genome and all other genomes in the given maf, using a tab delimited format.
Output table format has fields: querySpecies	targetSpecies	lengthOfQueryGenome	coverage	n-coverages (if specified)
For a pair of genomes A and B, the coverage of B on A is the proportion of sites in A that align to a base in B.
The n-coverage of B on A is the proportion of sites in A that align to n or more sites in B.
Options:
  -h, --help             show this help message and exit.
  -m, --maf              path to maf file.
  -s, --speciesOrChr     species or species.chromosome name, e.g. `hg19' or 'hg19.chr1',
                         if not specified reports results for every possible species.wildcard at
                         the end.
  -n, --nCoverage        report all n-coverages, for 1 <= n <= 128 instead of just
                         for n=1 (the default).
  -i, --identity         report coverage of identical bases.
  -l, --logLevel         Set logging level, either 'CRITICAL'/'INFO'/'DEBUG'.
  -a, --ignoreSpecies    Do all chromosomes-against-all-chromosomes coverage.
```


## Example
    $ mafCoverage --maf path/to/maf.maf --speciesOrChr hg19 ...
    ...

