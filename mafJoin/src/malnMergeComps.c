#include "malnMergeComps.h"
#include "malnBlkCursor.h"
#include "malnSet.h"
#include "malnBlk.h"
#include "malnComp.h"
#include "malnBlkSet.h"
#include "sonLibList.h"
#include <stdbool.h>

/* Determine if two components can be merged. That is, the sequences are
 * adjacent and are overlap or interleaved in the alignment */
static bool canMerge(struct malnComp *comp1, struct malnComp *comp2) {
    if (!malnComp_adjacentStrand(comp1, comp2)) {
        return false;
    } else {
        struct malnCompSeg *lastSeg1 = malnComp_getLastSeg(comp1);
        struct malnCompSeg *lastSeg2 = malnComp_getLastSeg(comp2);
        return ((malnCompSeg_getAlnEnd(lastSeg1) <= comp2->segs->alnStart)
                && (malnCompSeg_getEnd(lastSeg1) <= comp2->segs->start))
            || ((malnCompSeg_getAlnEnd(lastSeg2) <= comp1->segs->alnStart)
                && (malnCompSeg_getEnd(lastSeg2) <= comp1->segs->start));
    }
}

/* fast check to see if can any components be merged? */
static bool anyToMerge(struct malnBlk *blk) {
    for (struct malnComp *comp1 = blk->comps; comp1 != NULL; comp1 = comp1->next) {
        for (struct malnComp *comp2 = comp1->next; comp2 != NULL; comp2 = comp2->next) {
            if (canMerge(comp1, comp2)) {
                return true;
            }
        }
    }
    return false;
}

/* merge sequences in adjacent alignment components into comp1 */
static void mergeCompSegs(struct malnComp *comp1, struct malnComp *comp2) {
    struct malnCompSeg *lastSeg1 = malnComp_getLastSeg(comp1);
    for (struct malnCompSeg *seg2 = comp2->segs; seg2 != NULL; seg2 = seg2->next) {
        assert(malnCompSeg_getEnd(lastSeg1) == seg2->start);
        lastSeg1->next = malnCompSeg_constructClone(seg2);
        lastSeg1 = lastSeg1->next;
    }
}

/* Merge two adjacent components */
static void mergeAdjacentComps(struct malnBlk *blk, struct malnComp *comp1, struct malnComp *comp2) {
    assert(comp1->alnWidth == comp2->alnWidth);

    // make comp1 before comp2 in alignment
    if (comp1->segs->alnStart > comp2->segs->alnStart) {
        struct malnComp *compHold = comp1;
        comp1 = comp2;
        comp2 = compHold;
    }

    // copy contents
    comp1->start = min(comp1->start, comp2->start);
    comp1->end = max(comp1->end, comp2->end);
    comp1->chromStart = min(comp1->chromStart, comp2->chromStart);
    comp1->chromEnd = max(comp1->chromEnd, comp2->chromEnd);

    // copy comp2 segs to comp1
    mergeCompSegs(comp1, comp2);

    // remove from structures
    malnBlk_unlink(blk, comp2);
    malnComp_destruct(comp2);
}

/* one pass over components to merge a part. */
static bool mergeAdjacentCompsPass(struct malnBlk *blk) {
    for (struct malnComp *comp1 = blk->comps; comp1 != NULL; comp1 = comp1->next) {
        for (struct malnComp *comp2 = comp1->next; comp2 != NULL; comp2 = comp2->next) {
            if (canMerge(comp1, comp2)) {
                mergeAdjacentComps(blk, comp1, comp2);
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
 * This will replace the block in the set. */
static void mergeWithin(struct malnSet *malnSet, struct malnBlk *blk, struct malnBlkSet *newBlks) {
    // clone block so we are not modifying malnSet while updating
    struct malnBlk *newBlk = malnBlk_constructClone(blk);
    mergeAllAdjacentComps(newBlk);
    malnBlkSet_add(newBlks, newBlk);
    malnSet_markAsDeleted(malnSet, blk);
}

/* Merge overlapping and adjacent components within blocks. */
void malnMergeComps_merge(struct malnSet *malnSet) {
    // delay adding new blocks due to iterator
    struct malnBlkSet *newBlks = malnBlkSet_construct();
    struct malnBlkSetIterator *iter = malnSet_getBlocks(malnSet);
    struct malnBlk *blk;
    while ((blk = malnBlkSetIterator_getNext(iter)) != NULL) {
        // quick check first before overhead of merge
        if (anyToMerge(blk)) {
            mergeWithin(malnSet, blk, newBlks);
        }
    }
    malnBlkSetIterator_destruct(iter);
    malnSet_addBlks(malnSet, newBlks);
    malnBlkSet_destruct(newBlks);
    malnSet_deleteDying(malnSet);
}
