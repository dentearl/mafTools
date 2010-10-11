#ifndef malnBlkCursor_h
#define malnBlkCursor_h
#include "malnCompCursor.h"

/*
 * Cursor on all rows of an alignment block.  Optionally place a
 * reference component first.
 */
struct malnBlkCursor {
    struct malnBlk *blk;
    int numRows;
    int alnWidth;
    int alnIdx;
    struct malnCompCursor *rows;  // array of objects, not pointers
};


/* construct a new cursor, if refComp is not null, then force it to be first. If subsetComps
 * is not NULL, it is a NULL terminated list of a subset of components to include
 * in the order specified. */
struct malnBlkCursor *malnBlkCursor_construct(struct malnBlk *blk, struct malnComp *refComp, struct malnComp **subsetComps);

/* destructor */
void malnBlkCursor_destruct(struct malnBlkCursor *bc);

/* set the alignment cursors to the specify alignment column */
static inline void malnBlkCursor_setAlignCol(struct malnBlkCursor *bc, int alnIdx) {
    for (int i = 0; i < bc->numRows; i++) {
        malnCompCursor_setAlignCol(&(bc->rows[i]), alnIdx);
    }
    bc->alnIdx = alnIdx;
}

/* increment cursor, return false at the end */
static inline boolean malnBlkCursor_incr(struct malnBlkCursor *bc) {
    bc->alnIdx++;
    for (int i = 0; i < bc->numRows; i++) {
        malnCompCursor_incr(&(bc->rows[i]));
        assert(bc->rows[i].alnIdx == bc->alnIdx);
    }
    return bc->alnIdx < bc->alnWidth;
}

#endif
