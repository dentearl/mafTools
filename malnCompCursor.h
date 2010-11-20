#ifndef malnCompCursor_h
#define malnCompCursor_h
#include "malnComp.h"
#include "genome.h"
#include <stdbool.h>
#include <assert.h>

// FIXME: cursor should start out before the first position

/*
 * Column Cursor on a row of an alignment.
 */
struct malnCompCursor {
    struct malnComp *comp; 
    int alnIdx;             // alignment index, maybe set to before start or at end
    int pos;                // position, or next position if not aligned, or end
    char base;              // base at positon
    bool isAligned;         // is it aligned
};

/* make a new cursor at the start of a block */
static inline struct malnCompCursor malnCompCursor_make(struct malnComp *comp) {
    struct malnCompCursor cc;
    cc.comp = comp;
    cc.alnIdx = -1;
    cc.pos = comp->start;
    cc.base = '\0';
    cc.isAligned = false;
    return cc;
}

/* is the current position aligned */
static inline bool malnCompCursor_isAligned(struct malnCompCursor *cc) {
    return cc->isAligned;
}

/* is the current position at the end (after alignment) */
static inline bool malnCompCursor_atEnd(struct malnCompCursor *cc) {
    return cc->alnIdx == malnComp_getWidth(cc->comp);
}

/* is the current position at component end (no more aligned) */
static inline bool malnCompCursor_atCompEnd(struct malnCompCursor *cc) {
    assert(cc->pos <= cc->comp->chromEnd);
    return (cc->pos == cc->comp->chromEnd);
}

/* get the next position */
static inline int malnCompCursor_getNextPos(struct malnCompCursor *cc) {
    return cc->isAligned ? cc->pos+1 : cc->pos;
}

/* increment the alignment cursor, return false at the end. */
static inline bool malnCompCursor_incr(struct malnCompCursor *cc) {
    assert(cc->alnIdx <= malnComp_getWidth(cc->comp));
    if (cc->alnIdx == malnComp_getWidth(cc->comp)) {
        return false;
    }
    if (cc->isAligned) {
        cc->pos++;
    }
    cc->alnIdx++;
    if (cc->alnIdx == malnComp_getWidth(cc->comp)) {
        // reached end
        cc->pos = cc->comp->end;
        cc->isAligned = false;
        return false;
    } else {
        cc->base = malnComp_getCol(cc->comp, cc->alnIdx); 
        cc->isAligned = isBase(cc->base);
        return true;
    }
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

/* set the alignment cursor to the specify strand specific sequence position,
 * which maybe the end.*/
static inline void malnCompCursor_setSeqPos(struct malnCompCursor *cc, int pos) {
    assert((cc->comp->start <= pos) && (pos <= cc->comp->end));
    if (cc->pos > pos) {
        *cc = malnCompCursor_make(cc->comp);  // reset to start
    }
    while (!((cc->pos == pos) && malnCompCursor_isAligned(cc))) {
        if (!malnCompCursor_incr(cc)) {
            break;
        }
    }
    assert(cc->pos == pos);
    assert(malnCompCursor_atEnd(cc) || malnCompCursor_isAligned(cc));
}

/* advance the alignment cursor to the next aligned position, return false at the end  */
static inline bool malnCompCursor_advanceToNextPos(struct malnCompCursor *cc) {
    assert(cc->alnIdx <= malnComp_getWidth(cc->comp));
    do {
        if (!malnCompCursor_incr(cc)) {
            return false;
        }
    } while (!cc->isAligned);
    return true;
}

#endif
