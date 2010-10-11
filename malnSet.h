#ifndef malnSet_h
#define malnSet_h
#include "sonLibTypes.h"
struct Genomes;
struct Genome;
struct malnSet;
struct malnBlk;
struct malnComp;
struct Seq;

/* construct an empty malnSet  */
struct malnSet *malnSet_construct(struct Genomes *genomes);

/* Construct a malnSet from a MAF file. defaultBranchLength is used to
 * assign branch lengths when inferring trees from the MAF. */
struct malnSet *malnSet_constructFromMaf(struct Genomes *genomes, char *mafFileName, double defaultBranchLength, struct Genome *treelessRootGenome);

/* destructor */
void malnSet_destruct(struct malnSet *malnSet);

/* get associated genomes object  */
struct Genomes *malnSet_getGenomes(struct malnSet *malnSet);

/* add a block to a malnSet */
void malnSet_addBlk(struct malnSet *malnSet, struct malnBlk *blk);

/* remove a single component from malnSet range map */
void malnSet_removeComp(struct malnSet *malnSet, struct malnComp *comp);

/* remove a block from malnSet */
void malnSet_removeBlk(struct malnSet *malnSet, struct malnBlk *blk);

/* remove a block from malnSet and free the block */
void malnSet_deleteBlk(struct malnSet *malnSet, struct malnBlk *blk);

/* get iterator of the blocks. Don't remove or add blocks while in motion. */
stSortedSetIterator *malnSet_getBlocks(struct malnSet *malnSet);

stList *malnSet_getOverlappingPendingComps(struct malnSet *malnSet, struct Seq *seq, int chromStart, int chromEnd, unsigned treeLocFilter);

/* assert some sanity checks on a set */
void malnSet_assert(struct malnSet *malnSet);

/* assert done flag is set on all blocks */
void malnSet_assertDone(struct malnSet *malnSet);

/* clear done flag on all blocks */
void malnSet_clearDone(struct malnSet *malnSet);

/* write a malnSet to a MAF file  */
void malnSet_writeMaf(struct malnSet *malnSet, char *mafFileName);

/* return the number of blocks in the set */
int malnSet_getNumBlocks(struct malnSet *malnSet);

#endif
