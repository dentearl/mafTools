# mafToFastaStitcher

15 October 2012

## Author
[Dent Earl](https://github.com/dentearl/)

## Description
mafToFastaStitcherStrander is a program to take a multiple alignment format (MAF) file, some sequences, and then stitch together the alignment into a single multiple sequence fasta (MFA) file.

As an aside, the intended output is just dissimilar enough to what is created by multiz's maf2fasta that this tool is necessary but the output is similar enough as to be frustrating.

## Dependencies
* sonLib https://github.com/benedictpaten/sonLib/

## Installation
1. Download the package.
2. <code>cd</code> into the directory.
3. Type <code>make</code>.

## Use
<code>mafToFastaStitcher --maf alignment.maf --seqs seq.fa[,seq2.fa,...] --outMfa output.mfa --breakpointPenalty 10  [options] </code>

### Options
* <code>-h, --help</code>   show this help message and exit.
* <code>--maf</code>   input alignment maf file.
* <code>--seqs</code>   comma separated list of fasta sequences. each fasta may contain multiple entries. all sequences in the input alignment must be accounted for with an element in a fasta.
* <code>--outMfa</code>   multiple sequence fasta output file.
* <code>--breakpointPenalty</code>   number of <code>N</code> characters to insert into a sequence when a breakpoint is detected.
* <code>--interstitialSequence</code>   maximum length of interstitial sequence to be added (from a fasta) into the fasta before a breakpoint is declared and the <code>--breakpointPenalty</code> number of <code>N</code>'s is added instead.
* <code>--outMaf</code>    optional output to single block maf in addition to multiple sequence fasta output.

## Example
    $ mafToFastaStitcher --maf alignment.maf --seqs seq.fa,seq2.fa --breakpointPenalty 5 --outMfa output.mfa 

## Detailed input and output example

# Input maf
    ## maf version=1
    
    a score=0.0 status=test.input
    s ref.chr1   10 10 + 100 ACGTACGTAC
    s seq1.chr@   0 10 + 100 AAAAAAAAAA
    s seq2.chr&  10  5 + 100 -----CCCCC
    s seq6.chr1  10  5 + 100 -----GGGGG
    s seq7.chr20  0  5 + 100 AAAAA-----
    
    a score=0.0 status=test.input
    s ref.chr1   20 10 + 100 GTACGTACGT
    s seq2.chr!!  5  5 + 100 CCCCC-----
    s seq3.chr0  20  5 + 100 -----GGGGG
    s seq6.chr1  22  5 + 100 GGGGG-----
    
    a score=0.0 status=test.input
    s ref.chr1   30 10 + 100 ACGTACGTAC
    s seq4.chr1   0  5 - 100 GG-----GGG
    s seq5.chr2   0 10 + 100 CCCCCCCCCC
    s seq7.chr20 42  5 + 100 -----AAAAA
    
# Input sequences
Here in a single file, but could be broken across multiple files
    > ref.chr1 
    ggggggggggACGTACGTACGTACGTACGTACGTACGTACgg
    > seq1.chr@
    AAAAAAAAAAgg
    > seq2.chr&
    aaaaaaaaaaCCCCCaa
    > seq2.chr!!
    aaaaaCCCCCaa
    > seq3.chr0
    aaaaaaaaaaaaaaaaaaaGGGGGaa
    seq4.chr1
    aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
    aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaCCCCC
    > seq6.chr1
    aaaaaaaaaGGGGGaaaaaaaGGGGGaa
    > seq7.chr20
    AAAAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAAAAATT
    
# Expected multiple sequence fasta output
Note that when a sequence is only represented with a single chromosome, that chromosome will persist (as in ref.chr1) but when multiple chromosomes are present in the MAF that they are collapsed together.
    > ref.chr1
    ACGTACGTAC------------GTACGTACGT------------------
    -------------------ACGTACGTAC
    > seq1
    AAAAAAAAAA-----------------AAAAA------------------
    -----------------------------
    > seq2
    -----CCCCCNNNNN-------CCCCC-----------------------
    -----------------------------
    > seq6
    -----GGGGG-----aaaaaaaGGGGG-----------------------
    -----------------------------
    > seq7
    AAAAA---------------------------gggggggggggggggggg
    ggggggggggggggggggg-----AAAAA
    > seq3
    ---------------------------GGGGG------------------
    -----------------------------
    > seq4
    --------------------------------------------------
    -------------------GG-----GGG
    > seq5
    --------------------------------------------------
    -------------------CCCCCCCCCC


# optional maf output
where <code>--breakpointPenalty</code> is 5 (as seen in seq2) and <code>--interstitialSequence</code> is *at least* 17, as seen in seq7 (the long string of <code>g</code>'s is pulled in from the fasta).
    a score=0.0 status=test.expected
    s ref.chr1 10 30 + 100 ACGTACGTAC------------GTACGTACGT-------------------------------------ACGTACGTAC
    s seq1      0 15 +  15 AAAAAAAAAA-----------------AAAAA-----------------------------------------------
    s seq2      0 15 +  15 -----CCCCCNNNNN-------CCCCC----------------------------------------------------
    s seq6      0 17 +  17 -----GGGGG-----aaaaaaaGGGGG----------------------------------------------------
    s seq7      0 47 +  47 AAAAA---------------------------ggggggggggggggggggggggggggggggggggggg-----AAAAA
    s seq3      0  5 +   5 ---------------------------GGGGG-----------------------------------------------
    s seq4      0  5 -   5 ---------------------------------------------------------------------GG-----GGG
    s seq5      0 10 +  10 ---------------------------------------------------------------------CCCCCCCCCC
    
