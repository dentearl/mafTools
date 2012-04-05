# mafCoveragePickles

14 November 2011

## Authors
[Dent Earl](https://github.com/dentearl/)

## Description
mafCoveragePickles is a collection of scripts that perform simple
coverage based analyses on multiple alignment format (maf) files. 
All analyses use the pickle format, hense the name.

Fundamental to this tool is a naming convention for the name field.
Sequences should be named "`species.chrom`". This is important because
the tool `mafCoveragePickleGapAnalysis` performs an analysis by 
looking at all alignments of an entire genome to a specific sequence
(i.e. chromosome) of a second genome. The `chrom` aspect of the name
may contain periods, the `chrom` is defined as the string of 
characters in the name field that follow the first period.

## Dependencies
* Python 2.6 or 2.7
* [matplotlib](http://matplotlib.sourceforge.net/) if you intend to create coverage plots using `mafCoveragePicklePlotter`

## Installation
1. Download the package.
2. <code>cd</code> into the directory.
3. Type <code>make</code>.

## Workflow
Below is an example workflow. 
1. Run `mafCoveragePickleCreator` on a maf to decompose the maf into the pickle format.
2. Use the pickle format as input to `mafCoveragePickleSubsetExtractor` to pull out just the sequences of interest.
3. Use a subseted pickle as input to `mafCoveragePickleGapAnalysis` to perform a coverage and gap analysis of the alignment of a query genome onto a target sequence (specific chromosome).
4. Use a subseted pickle as input to `mafCoveragePicklePlotter` to create a dot plot between two sequences.

## Components
* `mafCoveragePickleCreator` - Takes as an input a maf file and creates a pickle format file. Output is pickle.
* `mafCoveragePickleGapAnalysis` - Takes as an input a pickle format file and a target and query sequence and performs an analysis to assess coverage and gap length distribution. Output is an xml.
* `mafCoveragePicklePlotter` - Takes as an input a pickle format file and a pair of sequences and produces a [dot-plot](http://en.wikipedia.org/wiki/Dot_plot_(bioinformatics)) showing the alignment between the two sequences. Output is pdf/eps/png.
* `mafCoveragePickleSubsetExtractor` - The pickle format files can be very large and can take a long time to read it can be beneficial to extract a subset of the pickle before performing analyses. This tool facilitates trimming the pickle. Output is pickle.

## Use
    mafCoveragePickleCreator.py --maf path/to/file.maf --species=species1,species2,... --pickle=output.pickle
    mafCoveragePickleGapAnalysis.py --pickle path/to/file.pickle --outfile=output.xml
    ...

## Test
<code>make test</code>
