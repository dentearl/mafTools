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

/* get reference genome object  */
struct Genome *malnSet_getRefGenome(struct malnSet *malnSet);

/* add a block to a malnSet */
void malnSet_addBlk(struct malnSet *malnSet, struct malnBlk *blk);

/* remove a block from malnSet */
void malnSet_removeBlk(struct malnSet *malnSet, struct malnBlk *blk);

/* remove a block from malnSet and free the block */
void malnSet_deleteBlk(struct malnSet *malnSet, struct malnBlk *blk);

/* get iterator of the blocks. Don't remove or add blocks while in motion. */
stSortedSetIterator *malnSet_getBlocks(struct malnSet *malnSet);

/* get list of slRefs to blks who's reference range overlaps the specified
 * range. */
stSortedSet *malnSet_getOverlapping(struct malnSet *malnSet, struct Seq *seq, int chromStart, int chromEnd);

/* clear done flag on all blocks */
void malnSet_clearDone(struct malnSet *malnSet);

/* write a malnSet to a MAF file  */
void malnSet_writeMaf(struct malnSet *malnSet, char *mafFileName);

#endif
