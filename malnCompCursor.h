#ifndef malnCompCursor_h
#define malnCompCursor_h
#include "malnComp.h"
#include "genome.h"
#include <stdbool.h>
#include <assert.h>

/*
 * Column Cursor on a row of an alignment.
 */
struct malnCompCursor {
    struct malnComp *comp; 
    int alnIdx;             // alignment index
    int pos;                // position, or previous if not aligned
    char base;              // base at positon
    bool isAligned;         // is it aligned
};

/* make a new cursor at the start of a block */
static inline struct malnCompCursor malnCompCursor_make(struct malnComp *comp) {
    struct malnCompCursor cc;
    cc.comp = comp;
    cc.alnIdx = 0;
    cc.pos = comp->start;
    cc.base = malnComp_getCol(comp, 0);
    cc.isAligned = isBase(cc.base);
    return cc;
}

/* increment the alignment cursor, return false at the end  */
static inline bool malnCompCursor_incr(struct malnCompCursor *cc) {
    assert(cc->alnIdx <= malnComp_getWidth(cc->comp));
    if (cc->alnIdx+1 >= malnComp_getWidth(cc->comp)) {
        if (cc->alnIdx < malnComp_getWidth(cc->comp)) {
            cc->alnIdx++;  // set to end
        }
        return false;
    }
    cc->alnIdx++;
    cc->base = malnComp_getCol(cc->comp, cc->alnIdx); 
    cc->isAligned = isBase(cc->base);
    if (cc->isAligned) {
        cc->pos++;
    }
    return true;
}

/* set the alignment cursor to the specify alignment column. Can be set to
 * the end of the alignment (to width) */
static inline void malnCompCursor_setAlignCol(struct malnCompCursor *cc, int alnIdx) {
    assert(alnIdx <= malnComp_getWidth(cc->comp));
    if (cc->alnIdx > alnIdx) {
        *cc = malnCompCursor_make(cc->comp);  // reset to start
    }
    while (cc->alnIdx < alnIdx) {
        malnCompCursor_incr(cc);
    }
}

/* set the alignment cursor to the specify strand specific sequence position */
static inline void malnCompCursor_setSeqPos(struct malnCompCursor *cc, int pos) {
    assert((cc->comp->start <= pos) && (pos < cc->comp->end));
    if (cc->pos > pos) {
        *cc = malnCompCursor_make(cc->comp);  // reset to start
    }
    while (cc->pos < pos) {
        malnCompCursor_incr(cc);
    }
}

/* advance the alignment cursor to the next aligned position, return false at the end  */
static inline bool malnCompCursor_nextPos(struct malnCompCursor *cc) {
    assert(cc->alnIdx < malnComp_getWidth(cc->comp));
    do {
        if (!malnCompCursor_incr(cc)) {
            return false;
        }
    } while (!cc->isAligned);
    return true;
}

/* is the current position aligned */
static inline bool malnCompCursor_isAligned(struct malnCompCursor *cc) {
    return cc->isAligned;
}

#endif
