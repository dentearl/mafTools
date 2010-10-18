#ifndef malnBlk_h
#define malnBlk_h
#include "mafJoinTypes.h"
#include <stdio.h>
#include <stdbool.h>
#include <stdarg.h>
#include <assert.h>
struct Seq;
struct Genome;

/* 
 * Multiple alignment block.
 */
struct malnBlk {
    int objId;                // unique, deterministic id.  Need because address space isn't deterministic
    struct malnSet *malnSet;  // set block is in, or NULL if not in a set.
    int alnWidth;
    struct malnComp *comps;   // components, root last
    mafTree *mTree;           // tree associated with block
    bool deleted;             // block is deleted, but not removed from structures yet.
};

/* constructor */
struct malnBlk *malnBlk_construct(void);

/* set the tree on a block */
void malnBlk_setTree(struct malnBlk *blk, mafTree *mTree);

/* constructor clone */
struct malnBlk *malnBlk_constructClone(struct malnBlk *srcBlk);

/* destructor */
void malnBlk_destruct(struct malnBlk *blk);

/* add a component */
void malnBlk_addComp(struct malnBlk *blk, struct malnComp *comp);

/* finish construction a block, setting component attributes and sorting
 * components */
void malnBlk_finish(struct malnBlk *blk);

/* Unlink a component from the block, also removing it from the tree and
 * malnSet map, if its in a set.  Component is not freed. */
void malnBlk_unlink(struct malnBlk *blk, struct malnComp *comp);

/* get the root component */
struct malnComp *malnBlk_getRootComp(struct malnBlk *blk);

/* find a component by seq and start, NULL if not found  */
struct malnComp *malnBlk_findBySeqStart(struct malnBlk *blk, struct Seq *seq, int start);

/* find a component by seq and chrom range, NULL if not found  */
struct malnComp *malnBlk_findCompByChromRange(struct malnBlk *blk, struct Seq *seq, int chromStart, int chromEnd);

/* block reverse complement */
struct malnBlk *malnBlk_reverseComplement(struct malnBlk *blk);

/* Shorten the range of a component in the block.  If the component is left
 * with no aligned bases, remove it.  Range must be at start or end of the
 * component. */
void malnBlk_shortenComp(struct malnBlk *blk, struct malnComp *comp, int newChromStart, int newChromEnd);

/* pad all components so they are the same width as the overall alignment. */
void malnBlk_pad(struct malnBlk *blk);

/* check for consistency within a components of a MAF, generating an
 * error if there are no consistent */
void malnBlk_validate(struct malnBlk *blk);

/* assert that the block is set-consistent */
void malnBlk_assert(struct malnBlk *blk);

/* compare two blocks for deterministic sorting */
int malnBlk_cmp(struct malnBlk *blk1, struct malnBlk *blk2);

/* If a block is in a set, mark it for deletion, otherwise just delete it */
void malnBlk_markOrDelete(struct malnBlk *blk);

/* construct an alignment block from a subrange of this block */
struct malnBlk *malnBlk_constrctSubrange(struct malnBlk *blk, int alnStart, int alnEnd);

/* print a block for debugging purposes */
void malnBlk_dumpv(struct malnBlk *blk, FILE *fh, const char *label, va_list args);

/* print a block for debugging purposes */
void malnBlk_dump(struct malnBlk *blk, FILE *fh, const char *label, ...)
#if defined(__GNUC__)
__attribute__((format(printf, 3, 4)))
#endif
;

#endif
