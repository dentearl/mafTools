#ifndef malnBlkCursor_h
#define malnBlkCursor_h
#include "malnCompCursor.h"

/*
 * Cursor on all rows of an alignment block.  Optionally place a
 * guide component first.
 */
struct malnBlkCursor {
    struct malnBlk *blk;
    int numRows;
    int alnWidth;
    int alnIdx;
    struct malnCompCursor *rows;  // array of objects, not pointers
};


/* construct a new cursor, if guideComp is not null, then force it to be first. If subsetComps
 * is not NULL, it is a NULL terminated list of a subset of components to include
 * in the order specified. */
struct malnBlkCursor *malnBlkCursor_construct(struct malnBlk *blk, struct malnComp *guideComp, struct malnComp **subsetComps);

/* destructor */
void malnBlkCursor_destruct(struct malnBlkCursor *bc);

/* set the alignment cursors to the specify alignment column */
static inline void malnBlkCursor_setAlignCol(struct malnBlkCursor *bc, int alnIdx) {
    for (int i = 0; i < bc->numRows; i++) {
        malnCompCursor_setAlignCol(&(bc->rows[i]), alnIdx);
    }
    bc->alnIdx = alnIdx;
}

/* is the current position at the end (after alignment) */
static inline bool malnBlkCursor_atEnd(struct malnBlkCursor *bc) {
    return bc->alnIdx == bc->alnWidth;
}

/* increment cursor, return false at the end.
 *  WARNING: cursor starts set at first position, not before */
static inline bool malnBlkCursor_incr(struct malnBlkCursor *bc) {
    bc->alnIdx++;
    for (int i = 0; i < bc->numRows; i++) {
        malnCompCursor_incr(&(bc->rows[i]));
        assert(bc->rows[i].alnIdx == bc->alnIdx);
    }
    return bc->alnIdx < bc->alnWidth;
}

#endif
