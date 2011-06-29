#mafJoin

**mafJoin** is a tool for combining pairs of [maf](http://genome.ucsc.edu/FAQ/FAQformat.html#format5) files that share a common sequence.

##Authors
Mark Diekhans

##Dependencies
* sonLib https://github.com/benedictpaten/sonLib/
* kentSource <code>git clone git://genome-source.cse.ucsc.edu/kent.git</code>

##Installation
1. Download the package. Consider making the parent of mafJoin a sibling directory to <code>sonLib</code>.
2. <code>cd</code> into the directory.
3. Type <code>make kentDir=/path/to/kent/src</code>.

##Use
<code>mafJoin [optional -treelessRoot1="sequence name" -treelessRoot2="sequence name" ...] "common sequence" first.maf second.maf out.maf</code>

Let there be two mafs, AC.maf and BC.maf, that share a common sequence C. A sequence in this context refers to the file wide sample name which is commonly in the format <code>species.chromosomeNumber</code>, e.g. sequence hg18 might have a line <code>hg18.chr19</code>. 

Let AC.maf contain two (in practice there may be any *n* >= 2) sequences, A and C and let BC.maf contain sequences B and C. A call to mafJoin would look like

<code>mafJoin -treelessRoot1=C -treelessRoot2=C C -maxBlkWidth=10000 -maxInputBlkWidth=1000 AC.maf BC.maf ABC.maf.tmp</code>

For the purposes in [evolverSimControl](https://github.com/dentearl/evolverSimControl/) the mafJoin command comes in two flavors. If we have a phylogeny (here in [Newick format](http://evolution.genetics.washington.edu/phylip/newicktree.html)) (((A, B)C, ) D); whereby A and B are siblings with parent node C and C is a child node of D then the way to create a single maf containing an alignment for ABCD is a two part process of first merging AC and BC into ABC and then merging ABC and CD into ABCD:

* <code>mafJoin -treelessRoot1=C -treelessRoot2=C C -maxBlkWidth=10000 -maxInputBlkWidth=1000 AC.maf BC.maf ABC.maf.tmp</code>
* <code>mv ABC.maf.tmp ABC.maf</code>
* <code>mafJoin -treelessRoot2=D C -maxBlkWidth=10000 -maxInputBlkWidth=1000 ABC.maf CD.maf ABCD.maf.tmp</code>

Note that in the second call to mafJoin we're only establishing the root (-treelessRoot) for the CD.maf sequence because mafJoin added tree information to the ABC.maf file making ABC.maf a tree maf (see docs/README for a technical defininition of a tree-maf). If you're using mafJoin in a multistep progressive merge then you only need to establish the -treelessRoot command if the maf in question lacks tree information. For example, imagine the tree (((A,B)C, (D, E)F)G; whereby G is the root node and has two children, C and F. C has two children A and B, and F has two children D and E. All of the joins would be:

* # ABC-G
* <code>mafJoin -treelessRoot1=C -treelessRoot2=C C -maxBlkWidth=10000 -maxInputBlkWidth=1000 AC.maf BC.maf ABC.maf.tmp</code>
* <code>mv ABC.maf.tmp ABC.maf</code>
* <code>mafJoin -treelessRoot2=G C -maxBlkWidth=10000 -maxInputBlkWidth=1000 ABC.maf CG.maf ABCG.maf.tmp</code>
* <code>mv ABCG.maf.tmp ABCG.maf</code>
* # DEF-G
* <code>mafJoin -treelessRoot1=F -treelessRoot2=F F -maxBlkWidth=10000 -maxInputBlkWidth=1000 DF.maf EF.maf DEF.maf.tmp</code>
* <code>mv DEF.maf.tmp DEF.maf</code>
* <code>mafJoin -treelessRoot2=D F -maxBlkWidth=10000 -maxInputBlkWidth=1000 DEF.maf FG.maf DEFG.maf.tmp</code>
* <code>mv DEFG.maf.tmp DEFG.maf</code>
* # ABCDEFG
* <code>mafJoin -maxBlkWidth=10000 -maxInputBlkWidth=1000 G ABCG.maf DEFG.maf ABCDEFG.maf.tmp</code>
* <code>mv ABCDEFG.maf.tmp ABCDEFG.maf</code>

Note that the final join does not have any -treelessRoot options since both incoming mafs are already tree-mafs.
