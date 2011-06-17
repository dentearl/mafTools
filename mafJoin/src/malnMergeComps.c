#include "malnMergeComps.h"
#include "malnBlkCursor.h"
#include "malnSet.h"
#include "malnBlk.h"
#include "malnComp.h"
#include "malnBlkSet.h"
#include "sonLibList.h"
#include <stdbool.h>

/* is a column in a state consistent with the start, or ignored because
 * neither component is aligned */
static bool isColConsistentState(struct malnCompCursor *cc0, struct malnCompCursor *cc1, bool aligned0, bool aligned1) {
    return ((!malnCompCursor_isAligned(cc0) && !malnCompCursor_isAligned(cc1))
            || ((malnCompCursor_isAligned(cc0) == aligned0) && (malnCompCursor_isAligned(cc1) == aligned1)));
}

/* find the next range of to components that is in at consistent state, that
 * is both aligned, only first aligned, or only second aligned.  Return
 * false when no more. Range is set to {-1,-1} when not aligned. */
static bool nextConsistentStateRange(struct malnBlkCursor *cursor, struct stRange *range0, struct stRange *range1) {
    *range0 = *range1 = stNullRange;
    if (malnBlkCursor_atEnd(cursor)) {
        return false;
    }
    struct malnCompCursor *cc0 = &(cursor->rows[0]);
    struct malnCompCursor *cc1 = &(cursor->rows[1]);
    bool aligned0 = malnCompCursor_isAligned(cc0);
    bool aligned1 = malnCompCursor_isAligned(cc1);
    if (aligned0) {
        range0->start = cc0->pos;
        range0->end = cc0->pos + 1;
    }
    if (aligned1) {
        range1->start = cc1->pos;
        range1->end = cc1->pos + 1;
    }
    
    // scan for consistent state, ignoring columns were neither are aligned
    while (malnBlkCursor_incr(cursor) && isColConsistentState(cc0, cc1, aligned0, aligned1)) {
        if (malnCompCursor_isAligned(cc0) || malnCompCursor_isAligned(cc1)) {
            if (aligned0) {
                range0->end = cc0->pos + 1;
            }
            if (aligned1) {
                range1->end = cc1->pos + 1;
            }
        }
    }
    return true;
}

/* check if two aligned component ranges are consistent, one maybe null */
static bool checkRangeConsistency(struct stRange prevAlignedRange0, struct stRange range0, struct stRange prevAlignedRange1, struct stRange range1) {
    if (!stRangeIsNull(range0) && !stRangeIsNull(range1)) {
        return stRangeEq(range0, range1);  // both aligned, must be the same
    }

    // ranges must consistently increase, otherwise they are interleaved
    if (!stRangeIsNull(range0) && !stRangeIsNull(prevAlignedRange1)) {
        if (!(range0.start >= prevAlignedRange1.end)) {
            return false;
        }
    }
    if (!stRangeIsNull(range1) && !stRangeIsNull(prevAlignedRange0)) {
        if (!(range1.start >= prevAlignedRange0.end)) {
            return false;
        }
    }
    return true;
}

/* scan two component in an alignment to determine if they are consistently
 * aligned and can be merged. */
static bool scanConsistentOverlapAlignment(struct malnBlk *blk, struct malnComp *comp1, struct malnComp *comp2) {
    struct malnComp *subsetComps[] = {comp1, comp2, NULL};
    struct malnBlkCursor *cursor = malnBlkCursor_construct(blk, NULL, subsetComps);

    struct stRange prevAlignedRange0 = stNullRange, prevAlignedRange1 = stNullRange;
    struct stRange range0 = stNullRange, range1 = stNullRange;
    bool isOk = true;
    while (isOk && nextConsistentStateRange(cursor, &range0, &range1)) {
        isOk = checkRangeConsistency(prevAlignedRange0, range0, prevAlignedRange1, range1);
        if (!stRangeIsNull(range0)) {
            prevAlignedRange0 = range0;
        }
        if (!stRangeIsNull(range1)) {
            prevAlignedRange1 = range1;
        }
    }
    malnBlkCursor_destruct(cursor);
    return isOk;
}

/* Check if two components are overlapping, adjacent and consistent in the
 * alignment.  That is, the same bases are aligned to the same columns and the
 * alignments are not interwoven */
static bool consistentOverlapAdjacent(struct malnBlk *blk, struct malnComp *comp1, struct malnComp *comp2) {
    if (!malnComp_overlapAdjacentStrand(comp1, comp2)) {
        return false;
    } else {
        return scanConsistentOverlapAlignment(blk, comp1, comp2);
    }
 }

/* upfront check for anything to merge */
static bool anyToMerge(struct malnBlk *blk) {
    for (struct malnComp *comp1 = blk->comps; comp1 != NULL; comp1 = comp1->next) {
        for (struct malnComp *comp2 = comp1->next; comp2 != NULL; comp2 = comp2->next) {
            if (consistentOverlapAdjacent(blk, comp1, comp2)) {
                return true;
            }
        }
    }
    return false;
}

/* bases in components being merged are consistent */
static void checkConsistentBases(struct malnComp *comp1, struct malnComp *comp2, int iCol) {
    if (isBase(malnComp_getCol(comp1, iCol)) && !baseEq(malnComp_getCol(comp1, iCol), malnComp_getCol(comp2, iCol))) {
        fprintf(stderr, "inconsistent sequences components being merged at column %d:\n", iCol);
        malnComp_dump(comp1, stderr, "comp1");
        malnComp_dump(comp2, stderr, "comp2");
        errAbort("invalid MAFs");
    }
}

/* merge sequences in overlapping alignment components into comp1 */
static void mergeCompSeqs(struct malnComp *comp1, struct malnComp *comp2) {
    for (int iCol = 0; iCol < malnComp_getWidth(comp2); iCol++) {
        if (isBase(malnComp_getCol(comp2, iCol))) {
            checkConsistentBases(comp1, comp2, iCol);
            malnComp_setCol(comp1, iCol, malnComp_getCol(comp2, iCol));
        }
    }
}

/* Merge two overlapping or adjacent components */
static void mergeOverlapAdjacentComps(struct malnBlk *blk, struct malnComp *comp1, struct malnComp *comp2) {
    assert(malnComp_getWidth(comp1) == malnComp_getWidth(comp2));

    // copy contents
    comp1->start = min(comp1->start, comp2->start);
    comp1->end = max(comp1->end, comp2->end);
    comp1->chromStart = min(comp1->chromStart, comp2->chromStart);
    comp1->chromEnd = max(comp1->chromEnd, comp2->chromEnd);
    mergeCompSeqs(comp1, comp2);
    if (malnComp_countAligned(comp1) != (comp1->end - comp1->start)) {  // FIXME
        fprintf(stderr, "comp1 aligned: %d, length %d\n", malnComp_countAligned(comp1), (comp1->end - comp1->start));
        malnComp_dump(comp2, stderr, "comp2");
        malnComp_dump(comp1, stderr, "comp1-new");
        malnComp_assert(comp1);
    }

    // remove from structures
    malnBlk_unlink(blk, comp2);
    malnComp_destruct(comp2);
}

/* one pass over components to merge a part. */
static bool mergeAdjacentCompsPass(struct malnBlk *blk) {
    for (struct malnComp *comp1 = blk->comps; comp1 != NULL; comp1 = comp1->next) {
        for (struct malnComp *comp2 = comp1->next; comp2 != NULL; comp2 = comp2->next) {
            if (consistentOverlapAdjacent(blk, comp1, comp2)) {
                mergeOverlapAdjacentComps(blk, comp1, comp2);
                return true;
            }
        }
    }
    return false;
}

/* Merge adjacent components.  Block should not be a member of a set,
 * since */
static void mergeAllAdjacentComps(struct malnBlk *blk) {
    // Repeated scan for pairs to merge until none are found, starting over
    // after a merge due to structure changing.
    bool foundPair;
    do {
        foundPair = mergeAdjacentCompsPass(blk);
    } while (foundPair);
    malnBlk_finish(blk);
}

/* Merge adjacent components within the specified block, if possible.
 * This will replace block in the set. */
static void mergeWithin(struct malnSet *malnSet, struct malnBlk *blk, struct malnBlkSet *newBlks) {
    if (anyToMerge(blk)) {
        // clone block so we are not modifying malnSet while updating
        struct malnBlk *mblk = malnBlk_constructClone(blk);
        mergeAllAdjacentComps(mblk);
        malnBlkSet_add(newBlks, mblk);
        malnSet_markAsDeleted(malnSet, blk);
    }
}

/* Merge adjacent components within the blocks of a set. */
void malnMergeComps_merge(struct malnSet *malnSet) {
    // delay add and delete due to iterator
    struct malnBlkSet *newBlks = malnBlkSet_construct();
    struct malnBlkSetIterator *iter = malnSet_getBlocks(malnSet);
    struct malnBlk *blk;
    while ((blk = malnBlkSetIterator_getNext(iter)) != NULL) {
        mergeWithin(malnSet, blk, newBlks);
    }
    malnBlkSetIterator_destruct(iter);
    malnSet_addBlks(malnSet, newBlks);
    malnBlkSet_destruct(newBlks);
    malnSet_deleteDying(malnSet);
}
