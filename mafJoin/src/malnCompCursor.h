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
    int alnWidth;           // width of the alignment
    int alnIdx;             // alignment index, maybe set to before start or at end
    int pos;                // position, or next position if not aligned, or end
    bool isAligned;         // is it aligned
    struct malnCompSeg *seg;  // current or next segment of component
};
// FIXME: is alnWidth really needed?

/* make a new cursor befor the start of a block */
static inline struct malnCompCursor malnCompCursor_make(struct malnComp *comp) {
    struct malnCompCursor cc;
    cc.comp = comp;
    cc.alnIdx = -1;
    cc.alnWidth = comp->alnWidth;
    cc.pos = comp->start;
    cc.isAligned = false;
    cc.seg = comp->segs;
    return cc;
}

/* assert that the cursor appears valid */
static inline void malnCompCursor_assert(struct malnCompCursor *cc) {
    assert(cc->alnWidth == cc->comp->alnWidth);
    assert((-1 <= cc->alnIdx) && (cc->alnIdx <= cc->comp->alnWidth));
    assert((cc->comp->start <= cc->pos) && (cc->pos <= cc->comp->end));
}


/* get the current column */
static inline char malnCompCursor_getCol(struct malnCompCursor *cc) {
    if (cc->isAligned) {
        return  malnCompSeg_getBase(cc->seg, cc->alnIdx);
    } else {
        return '-';
    }
}

/* is the current position aligned */
static inline bool malnCompCursor_isAligned(struct malnCompCursor *cc) {
    return cc->isAligned;
}

/* is the current position at the end (after alignment) */
static inline bool malnCompCursor_atEnd(struct malnCompCursor *cc) {
    return cc->alnIdx == cc->alnWidth;
}

/* is the current position at component end (no more aligned) */
static inline bool malnCompCursor_atCompEnd(struct malnCompCursor *cc) {
    malnCompCursor_assert(cc);
    return (cc->pos == cc->comp->end);
}

/* get the next position */
static inline int malnCompCursor_getNextPos(struct malnCompCursor *cc) {
    return cc->isAligned ? cc->pos+1 : cc->pos;
}

/* advance cursor to the next seg.  If advance past the end, set past the
 * end of the last block. */
static inline bool malnCompCursor_nextSeg(struct malnCompCursor *cc) {
    malnCompCursor_assert(cc);
    if (cc->seg == NULL) {
        return false;
    } else if (cc->seg->next == NULL) {
        cc->alnIdx = malnCompSeg_getAlnEnd(cc->seg);
        cc->pos = cc->comp->end;
        cc->isAligned = false;
        cc->seg = NULL;
        return false;
    } else {
        cc->seg = cc->seg->next;
        cc->alnIdx = cc->seg->alnStart;
        cc->pos = cc->seg->start;
        cc->isAligned = true;
        return true;
    }
}

/* increment the alignment cursor, return false at the end. */
static inline bool malnCompCursor_incr(struct malnCompCursor *cc) {
    malnCompCursor_assert(cc);
    if (cc->alnIdx == cc->alnWidth) {
        // already at the end
        return false;
    }
    if (cc->isAligned) {
        cc->pos++;
    }
    cc->alnIdx++;
    if (cc->alnIdx == cc->alnWidth) {
        // reached end of alignment
        cc->pos = cc->comp->end;
        cc->isAligned = false;
        cc->seg = NULL;
        return false;
    } else {
        if ((cc->seg != NULL) && (cc->alnIdx == malnCompSeg_getAlnEnd(cc->seg))) {
            // advance to next seg
            cc->seg = cc->seg->next;
        }
        cc->isAligned = (cc->seg != NULL) && malnCompSeg_alnIdxInSeg(cc->seg, cc->alnIdx);
        return true;
    }
}

/* Set the cursor to an alignment index within the current seg */
static inline void malnCompCursor_setAlnIdxForSeg(struct malnCompCursor *cc, int alnIdx) {
    assert((cc->seg != NULL) && malnCompSeg_alnIdxInSeg(cc->seg, alnIdx));
    cc->alnIdx = alnIdx;
    cc->pos = cc->seg->start + (alnIdx - cc->seg->alnStart);
    cc->isAligned = true;
}


/* set the alignment cursor to the specify alignment column. Can be set to
 * the end of the alignment (to width) */
static inline void malnCompCursor_setAlignCol(struct malnCompCursor *cc, int alnIdx) {
    assert(alnIdx <= cc->comp->alnWidth);
    if (cc->alnIdx > alnIdx) {
        *cc = malnCompCursor_make(cc->comp);  // reset to start
    }
    while (cc->seg != NULL) {
        if (alnIdx < cc->seg->alnStart) {
            cc->alnIdx = alnIdx;
            cc->isAligned = false;
            break;
        } else if (malnCompSeg_alnIdxInSeg(cc->seg, alnIdx)) {
            malnCompCursor_setAlnIdxForSeg(cc, alnIdx);
            break;
        } else {
            if (!malnCompCursor_nextSeg(cc)) {
                break;
            }
        }
    }
    if (cc->seg == NULL) {
        // after last seg
        cc->alnIdx = alnIdx;
        cc->isAligned = false;
    }
}

/* Set the cursor to an position within the current seg */
static inline void malnCompCursor_setPosForSeg(struct malnCompCursor *cc, int pos) {
    assert((cc->seg != NULL) && malnCompSeg_containsPos(cc->seg, pos));
    cc->pos = pos;
    cc->alnIdx = cc->seg->alnStart + (pos - cc->seg->start);
    cc->isAligned = true;
}

/* Set the cursor to position at the end */
static inline void malnCompCursor_setToEnd(struct malnCompCursor *cc) {
    cc->pos = cc->comp->end;
    cc->alnIdx = cc->alnWidth;
    cc->isAligned = false;
    cc->seg = NULL;
}

/* set the alignment cursor to the specify strand specific sequence position,
 * which maybe the end.*/
static inline void malnCompCursor_setSeqPos(struct malnCompCursor *cc, int pos) {
    assert(malnComp_containsPos(cc->comp, pos) || (pos == cc->comp->end));
    if (pos == cc->comp->end) {
        malnCompCursor_setToEnd(cc);
    } else {
        if (cc->pos > pos) {
            *cc = malnCompCursor_make(cc->comp);  // reset to start
        }
        while ((cc->seg != NULL) && !malnCompSeg_containsPos(cc->seg, pos)) {
            malnCompCursor_nextSeg(cc);
        }
        malnCompCursor_setPosForSeg(cc, pos);
        assert((cc->pos == pos) && malnCompCursor_isAligned(cc));
    }
}

/* If the cursor is not aligned, advance it to next aligned position, return false at the end  */
static inline bool malnCompCursor_advanceToAligned(struct malnCompCursor *cc) {
    malnCompCursor_assert(cc);
    if (!cc->isAligned) {
        if (cc->seg == NULL) {
            return false;
        }
        cc->alnIdx = cc->seg->alnStart;
        cc->pos = cc->seg->start;
        cc->isAligned = true;
    }
    return true;
}

/* print cursor info for debugging */
static inline void malnCompCursor_dump(struct malnCompCursor *cc, FILE *fh) {
    fprintf(fh, "idx: %d pos %d %s\n", cc->alnIdx, cc->pos, (cc->isAligned ? "Y" : "N"));
}
#endif
