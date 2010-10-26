#ifndef malnSet_h
#define malnSet_h
#include <stdarg.h>
#include "sonLibTypes.h"
struct Genomes;
struct Genome;
struct malnSet;
struct malnBlk;
struct malnComp;
struct Seq;
struct malnBlkSet;

/* construct an empty malnSet  */
struct malnSet *malnSet_construct(struct Genomes *genomes);

/* Construct a malnSet from a MAF file. defaultBranchLength is used to
 * assign branch lengths when inferring trees from the MAF. */
struct malnSet *malnSet_constructFromMaf(struct Genomes *genomes, char *mafFileName, int maxInputBlkWidth, double defaultBranchLength, struct Genome *treelessRootGenome);

/* destructor */
void malnSet_destruct(struct malnSet *malnSet);

/* get associated genomes object  */
struct Genomes *malnSet_getGenomes(struct malnSet *malnSet);

/* add a block to a malnSet */
void malnSet_addBlk(struct malnSet *malnSet, struct malnBlk *blk);

/* add all blocks in to a malnSet */
void malnSet_addBlks(struct malnSet *malnSet, struct malnBlkSet *blks);

/* add a component to the range map. */
void malnSet_addComp(struct malnSet *malnSet, struct malnComp *comp);

/* remove a single component from malnSet range map */
void malnSet_removeComp(struct malnSet *malnSet, struct malnComp *comp);

/* remove a block from malnSet */
void malnSet_removeBlk(struct malnSet *malnSet, struct malnBlk *blk);

/* pop a block from the set.  order block are returned is arbitrary */
struct malnBlk *malnSet_popBlk(struct malnSet *malnSet);

/* get iterator of the blocks. Don't remove or add blocks while in motion. */
struct malnBlkSetIterator *malnSet_getBlocks(struct malnSet *malnSet);

/* Get a table with the set of blocks. Useful when one needs to add remove or
 * add blocks while scanning. */
struct malnBlkSet *malnSet_getBlockSetCopy(struct malnSet *malnSet);

/* Get a list of components that overlap the specified guide range and are
 * in blocks not flagged as done or dying, passing treeLoc filters, and not in
 * option doneBlks.  Return NULL if no overlaps. */
stList *malnSet_getOverlappingPendingComps(struct malnSet *malnSet, struct Seq *seq, int chromStart, int chromEnd, unsigned treeLocFilter, struct malnBlkSet *doneBlks);

/* Get a list of components that overlap or are adjacent to the specified
 * guide range and are in blocks not flagged as done or dying, passing
 * treeLoc filters, and not in option doneBlks.  Return NULL if no
 * overlaps. */
stList *malnSet_getOverlappingAdjacentPendingComps(struct malnSet *malnSet, struct Seq *seq, int chromStart, int chromEnd, unsigned treeLocFilter, struct malnBlkSet *doneBlks);

/* validate set for consistency */
void malnSet_validate(struct malnSet *malnSet);

/* assert some sanity checks on a set */
void malnSet_assert(struct malnSet *malnSet);

/* Record a block as deleted.  It is allowed to add blocks with are not
 * associated with a malnSet, as a way of cleaning up intermediate blocks. */
void malnSet_markAsDeleted(struct malnSet *malnSet, struct malnBlk *blk);

/* delete blocks marked as dying */
void malnSet_deleteDying(struct malnSet *malnSet);

/* write a malnSet to a MAF file  */
void malnSet_writeMaf(struct malnSet *malnSet, char *mafFileName);

/* print set for debugging */
void malnSet_dumpv(struct malnSet *malnSet, FILE *fh, const char *label, va_list args);

/* print set for debugging */
void malnSet_dump(struct malnSet *malnSet, FILE *fh, const char *label, ...)
#if defined(__GNUC__)
__attribute__((format(printf, 3, 4)))
#endif
;

#endif
