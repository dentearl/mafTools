#ifndef stMalnCompCursor_h
#define stMalnCompCursor_h
#include "stMalnComp.h"
#include "genome.h"
#include <stdbool.h>
#include <assert.h>

/*
 * Column Cursor on a row of an alignment.
 */
struct stMalnCompCursor {
    struct stMalnComp *comp; 
    int alnIdx;             // alignment index
    int pos;                // position, or previous if not aligned
    char base;              // base at positon
    bool isAligned;         // is it aligned
};

/* make a new cursor at the start of a block */
static inline struct stMalnCompCursor stMalnCompCursor_make(struct stMalnComp *comp) {
    struct stMalnCompCursor cc;
    cc.comp = comp;
    cc.alnIdx = 0;
    cc.pos = comp->start;
    cc.base = stMalnComp_getCol(comp, 0);
    cc.isAligned = isBase(cc.base);
    return cc;
}

/* increment the alignment cursor, return false at the end  */
static inline bool stMalnCompCursor_incr(struct stMalnCompCursor *cc) {
    assert(cc->alnIdx <= stMalnComp_getWidth(cc->comp));
    if (cc->alnIdx+1 >= stMalnComp_getWidth(cc->comp)) {
        if (cc->alnIdx < stMalnComp_getWidth(cc->comp)) {
            cc->alnIdx++;  // set to end
        }
        return false;
    }
    cc->alnIdx++;
    cc->base = stMalnComp_getCol(cc->comp, cc->alnIdx); 
    cc->isAligned = isBase(cc->base);
    if (cc->isAligned) {
        cc->pos++;
    }
    return true;
}

/* set the alignment cursor to the specify alignment column. Can be set to
 * the end of the alignment (to width) */
static inline void stMalnCompCursor_setAlignCol(struct stMalnCompCursor *cc, int alnIdx) {
    assert(alnIdx <= stMalnComp_getWidth(cc->comp));
    if (cc->alnIdx > alnIdx) {
        *cc = stMalnCompCursor_make(cc->comp);  // reset to start
    }
    while (cc->alnIdx < alnIdx) {
        stMalnCompCursor_incr(cc);
    }
}

/* set the alignment cursor to the specify strand specific sequence position */
static inline void stMalnCompCursor_setSeqPos(struct stMalnCompCursor *cc, int pos) {
    assert((cc->comp->start <= pos) && (pos < cc->comp->end));
    if (cc->pos > pos) {
        *cc = stMalnCompCursor_make(cc->comp);  // reset to start
    }
    while (cc->pos < pos) {
        stMalnCompCursor_incr(cc);
    }
}

/* advance the alignment cursor to the next aligned position, return false at the end  */
static inline bool stMalnCompCursor_nextPos(struct stMalnCompCursor *cc) {
    assert(cc->alnIdx < stMalnComp_getWidth(cc->comp));
    do {
        if (!stMalnCompCursor_incr(cc)) {
            return false;
        }
    } while (!cc->isAligned);
    return true;
}

/* is the current position aligned */
static inline bool stMalnCompCursor_isAligned(struct stMalnCompCursor *cc) {
    return cc->isAligned;
}

#endif
