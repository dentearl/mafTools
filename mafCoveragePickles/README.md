#mafIndelDistribution.py

14 November 2011

##Authors
[Dent Earl](https://github.com/dentearl/)

##Description
mafIndelDistribution is a script that operates on a single maf
(multiple alignment format) file and extracts information from
the alignments of pairs of sequences. Output is xml which will
contain a 'gaps' tag that contains a comma separated list of all
indels within the maf for specified pairs. Additionally, 
information on inter-chromosome coverage between pairs is 
provided. To this end, the maf sequence name field must have the
format
species.chromosome
e.g.
hg19.chr22

For slightly more detailed examples, check the 
test.mafIndelDistribution.py unittest.

##Dependencies
* Python 2.6 &le; version &lt; 3.0

##Installation
1. Download the package.
2. <code>cd</code> into the directory.
3. Type <code>make</code>.

##Use
<code>mafIndelDistribution.py --maf path/to/file.maf --species=species1,species2,...</code>

###Options
* <code>mafIndelDistribution.py</code>
* <code>-h, --help</code>         show this help message and exit
* <code>--maf=MAF</code>          input maf file
* <code>--species=SPECIES</code>  comma separated list of species names to include in output
* <code>--outfile=OUTFILE</code>  directory where outfile will be written default = summary.xml

##Test
<code>make test</code>
