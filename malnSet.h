#ifndef malnSet_h
#define malnSet_h
#include "sonLibTypes.h"
struct Genomes;
struct Genome;
struct malnSet;
struct malnBlk;
struct Seq;

/* construct an empty malnSet  */
struct malnSet *malnSet_construct(struct Genomes *genomes, struct Genome *refGenome);

/* Construct a malnSet from a MAF file. defaultBranchLength is used to
 * assign branch lengths when inferring trees from pair-wise MAFs. */
struct malnSet *malnSet_constructFromMaf(struct Genomes *genomes, struct Genome *refGenome, char *mafFileName, double defaultBranchLength);

/* get associated genomes object  */
struct Genomes *malnSet_getGenomes(struct malnSet *malnSet);

/* add a block to a malnSet */
void malnSet_addBlk(struct malnSet *malnSet, struct malnBlk *blk);

/* write a malnSet to a MAF file  */
void malnSet_writeMaf(struct malnSet *malnSet, char *mafFileName);

/* get list of blocks.  Don't modify this list. */
struct malnBlk *malnSet_getBlocks(struct malnSet *malnSet);

/* get list of slRefs to blks who's reference range overlaps the specified
 * range.  slFreeList of results when done. */
struct slRef *malnSet_getOverlapping(struct malnSet *malnSet, struct Seq *seq, int chromStart, int chromEnd);

#endif
