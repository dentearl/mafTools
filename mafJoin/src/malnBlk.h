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
    int alnWidth;             // during construction, this will have the max width of any component
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

/* Clear the sequences for a dying block.  This supports freeing up large
 * amounts of memory immediately without having to update the metadata,
 * requiring iterators to be invalidated, etc. */
void malnBlk_freeSeqMem(struct malnBlk *blk);

/* If a block is part of a set, mark it for deletion in that set,
 * otherwise deleted immediately*/
void malnBlk_markOrDelete(struct malnBlk *blk);

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

/* pad all components so they are the same width as the overall alignment. */
void malnBlk_pad(struct malnBlk *blk);

/* check for consistency within a components of a MAF, generating an
 * error if there are no consistent */
void malnBlk_validate(struct malnBlk *blk);

/* assert that the block is set-consistent.  ncLinkCheck can be false for
 * sanity check before links and tree are constructed. */
void malnBlk_assert(struct malnBlk *blk, bool ncLinkCheck);

/* compare two blocks for deterministic sorting */
int malnBlk_cmp(struct malnBlk *blk1, struct malnBlk *blk2);

/* construct an alignment block from a subrange of this block */
struct malnBlk *malnBlk_constructSubrange(struct malnBlk *blk, int alnStart, int alnEnd);

enum {
    malnBlk_dumpDefault   = 0x00,
    malnBlk_dumpNoSeqs    = 0x01,
    malnBlk_dumpInclTree  = 0x02,
};

/* print a block for debugging purposes */
void malnBlk_dumpv(struct malnBlk *blk, FILE *fh, unsigned opts, const char *label, va_list args);

/* print a block for debugging purposes */
void malnBlk_dump(struct malnBlk *blk, FILE *fh, unsigned opts, const char *label, ...)
#if defined(__GNUC__)
__attribute__((format(printf, 4, 5)))
#endif
;

#endif
