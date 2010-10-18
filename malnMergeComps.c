#include "malnMergeComps.h"
#include "malnBlkCursor.h"
#include "malnSet.h"
#include "malnBlk.h"
#include "malnComp.h"
#include "malnBlkSet.h"
#include "sonLibList.h"
#include <stdbool.h>

static bool debug = false;  // FIXME: tmp

/* FIXME: it would be more efficient to generalize malnJoinBlks to handle this.
 */

/* fast check for anything to merge */
static bool anyToMerge(struct malnBlk *blk) {
    for (struct malnComp *comp1 = blk->comps; comp1 != NULL; comp1 = comp1->next) {
        for (struct malnComp *comp2 = comp1->next; comp2 != NULL; comp2 = comp2->next) {
            if (malnComp_overlapAdjacentStrand(comp1, comp2)) {
                return true;
            }
        }
    }
    return false;
}

/* check for column consistency */
static inline bool checkColConsistency(struct malnBlkCursor *cursor, int *prevPos) {
    struct malnCompCursor *cc0 = &(cursor->rows[0]);
    struct malnCompCursor *cc1 = &(cursor->rows[1]);
    bool aligned0 = malnCompCursor_isAligned(cc0);
    bool aligned1 = malnCompCursor_isAligned(cc1);
    if (!(aligned0 || aligned1)) {
        return true; // neither aligned
    }
    if (aligned0 && aligned1) {
        *prevPos = cc0->pos;
        return (cc0->pos == cc1->pos);  // must be same position
    }
    struct malnCompCursor *alnCc = aligned0 ? cc0 : cc1;
    if ((*prevPos >= 0) && (alnCc->pos < *prevPos)) {
        return false; // interleaved
    }
    *prevPos = alnCc->pos;
    return true;  // one aligned
}

/* Check if two components are overlapping, adjacent and consistent in the
 * alignment.  That is, the same bases are aligned to the same columns and the
 * alignments are not interwoven */
static bool consistentOverlapAdjacent(struct malnBlk *blk, struct malnComp *comp1, struct malnComp *comp2) {
    if (!malnComp_overlapAdjacentStrand(comp1, comp2)) {
        return false;
    }
    struct malnComp *subsetComps[] = {comp1, comp2, NULL};
    struct malnBlkCursor *cursor = malnBlkCursor_construct(blk, NULL, subsetComps);
    bool isOk = true;
    int prevPos = -1;
    while (isOk && malnBlkCursor_incr(cursor)) {
        isOk = checkColConsistency(cursor, &prevPos);
    }
    malnBlkCursor_destruct(cursor);
    return isOk;
}

/* merge sequences in overlapping alignment components into comp1 */
static void mergeCompSeqs(struct malnComp *comp1, struct malnComp *comp2) {
    for (int i = 0; i < malnComp_getWidth(comp2); i++) {
        if (isBase(malnComp_getCol(comp2, i))) {
            assert((!isBase(malnComp_getCol(comp1, i))) || baseEq(malnComp_getCol(comp1, i), malnComp_getCol(comp2, i)));
            malnComp_setCol(comp1, i, malnComp_getCol(comp2, i));
        }
    }
}

/* Merge two adjacent components */
static void mergeAdjacentComps(struct malnBlk *blk, struct malnComp *comp1, struct malnComp *comp2) {
    assert(malnComp_getWidth(comp1) == malnComp_getWidth(comp2));

    // copy contents
    comp1->start = min(comp1->start, comp2->start);
    comp1->end = max(comp1->end, comp2->end);
    comp1->chromStart = min(comp1->chromStart, comp2->chromStart);
    comp1->chromEnd = max(comp1->chromEnd, comp2->chromEnd);
    mergeCompSeqs(comp1, comp2);

    // remove from structures
    malnBlk_unlink(blk, comp2);
    malnComp_destruct(comp2);
}

/* one pass over components to merge a part */
static bool mergeAdjacentCompsPass(struct malnBlk *blk) {
    malnBlk_assert(blk); // FIXME: tmp
    for (struct malnComp *comp1 = blk->comps; comp1 != NULL; comp1 = comp1->next) {
        for (struct malnComp *comp2 = comp1->next; comp2 != NULL; comp2 = comp2->next) {
            if (consistentOverlapAdjacent(blk, comp1, comp2)) {
                mergeAdjacentComps(blk, comp1, comp2);
                malnBlk_assert(blk); // FIXME: tmp
                return true;
            }
        }
    }
    return false;
}

/* Merge adjacent components.  Block should not be a member of a set,
 * since */
static void mergeAllAdjacentComps(struct malnBlk *blk) {
    assert(blk->malnSet == NULL);

    // repeated scan for pairs to merge until none are found, starting over
    // after a merge due to structure changing
    bool foundPair;
    do {
        foundPair = mergeAdjacentCompsPass(blk);
    } while (foundPair);
}

/* Merge adjacent components within the specified block, if possible.
 * This will replace block in the set. */
static void mergeWithin(struct malnSet *malnSet, struct malnBlk *blk, struct malnBlkSet *newBlks) {
    if (anyToMerge(blk)) {
        struct malnBlk *mblk = malnBlk_constructClone(blk);
        malnBlk_assert(mblk); // FIXME: tmp
        if (debug) {
            malnBlk_dump(mblk, stderr, "mergeWithin:before");
        }
        mergeAllAdjacentComps(mblk);
        malnBlkSet_add(newBlks, mblk);
        malnBlk_markOrDelete(blk);
        if (debug) {
            malnBlk_dump(mblk, stderr, "mergeWithin:after");
        }
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
