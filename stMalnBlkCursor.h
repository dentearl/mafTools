#ifndef stMalnBlkCursor_h
#define stMalnBlkCursor_h
#include "stMalnCompCursor.h"

/*
 * Cursor on all rows of an alignment block.  Optionally place a
 * reference component first.
 */
struct stMalnBlkCursor {
    struct stMalnBlk *blk;
    int numRows;
    int alnWidth;
    int alnIdx;
    struct stMalnCompCursor *rows;  // array of objects, not pointers
};


/* construct a new cursor, if refComp is not null, then force it to be first */
struct stMalnBlkCursor *stMalnBlkCursor_construct(struct stMalnBlk *blk, struct stMalnComp *refComp);

/* destructor */
void stMalnBlkCursor_destruct(struct stMalnBlkCursor *bc);

/* set the alignment cursors to the specify alignment column */
static inline void stMalnBlkCursor_setAlignCol(struct stMalnBlkCursor *bc, int alnIdx) {
    for (int i = 0; i < bc->numRows; i++) {
        stMalnCompCursor_setAlignCol(&(bc->rows[i]), alnIdx);
    }
    bc->alnIdx = alnIdx;
}

/* increment cursor, return false at the end */
static inline boolean stMalnBlkCursor_incr(struct stMalnBlkCursor *bc) {
    bc->alnIdx++;
    for (int i = 0; i < bc->numRows; i++) {
        stMalnCompCursor_incr(&(bc->rows[i]));
        assert(bc->rows[i].alnIdx == bc->alnIdx);
    }
    return bc->alnIdx < bc->alnWidth;
}

#endif
