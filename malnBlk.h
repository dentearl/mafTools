#ifndef malnBlk_h
#define malnBlk_h
#include "mafJoinTypes.h"
#include <stdio.h>
#include <stdbool.h>
struct Seq;
struct Genome;

/* 
 * Multiple alignment block.
 */
struct malnBlk {
    struct malnSet *malnSet;  // set block is in, or NULL if not in a set.
    int alnWidth;
    struct malnComp *comps;   // components
    mafTree *mTree;           // tree associated with block
    bool done;                // finished some kind of processing
};

/* constructor */
struct malnBlk *malnBlk_construct(mafTree *mTree);

/* constructor clone */
struct malnBlk *malnBlk_constructClone(struct malnBlk *srcBlk);

/* destructor */
void malnBlk_destruct(struct malnBlk *blk);

/* add a component */
void malnBlk_addComp(struct malnBlk *blk, struct malnComp *comp);

/* finish construction a block, setting component attributes and sorting
 * components */
void malnBlk_finish(struct malnBlk *blk);

/* get the root component */
struct malnComp *malnBlk_getRootComp(struct malnBlk *blk);

/* find a component by seq and start, NULL if not found  */
struct malnComp *malnBlk_findBySeqStart(struct malnBlk *blk, struct Seq *seq, int start);

/* find a component by seq and chrom range, NULL if not found  */
struct malnComp *malnBlk_findCompByChromRange(struct malnBlk *blk, struct Seq *seq, int chromStart, int chromEnd);

/* block reverse complement */
struct malnBlk *malnBlk_reverseComplement(struct malnBlk *blk);

/* pad all components so they are the same width as the overall alignment. */
void malnBlk_pad(struct malnBlk *blk);

/* assert that the block is set-consistent */
void malnBlk_assert(struct malnBlk *blk);

/* print a block for debugging purposes */
void malnBlk_dump(struct malnBlk *blk, const char *label, FILE *fh);

#endif
