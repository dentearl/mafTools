#include "malnJoinWithinSet.h"
#include "malnSet.h"
#include "malnBlk.h"
#include "malnComp.h"
#include "malnBlkSet.h"
#include "malnJoinBlks.h"
#include "sonLibSortedSet.h"
#include "sonLibList.h"
#include "genome.h"
#include "common.h"
#include <stdbool.h>
#include <unistd.h>

/* join two blocks by associated root components, create a third
 * block. Return that new block. Joined blocks marked as deleted */
static struct malnBlk *joinBlksAtRoot(struct malnSet *malnSet, struct malnBlk *blk1, struct malnBlk *blk2) {
    struct malnBlk *joinBlk = malnJoinBlks(malnBlk_getRootComp(blk1), malnBlk_getRootComp(blk2));
    malnSet_markAsDeleted(malnSet, blk1);
    malnSet_markAsDeleted(malnSet, blk2);
    return joinBlk;
}

/* Attempt to join a block with other blocks. Return updated block when one
 * join is achieved, or NULL if no join was done. */
static struct malnBlk *joinBlkWithDupsPass(struct malnSet *malnSet, struct malnBlk *joinBlk) {
    // n.b. don't include adjacent, as that can causes very large blocks at
    // the root.
    struct malnBlk *newJoinBlk =  NULL;
    struct malnComp *joinComp = malnBlk_getRootComp(joinBlk);
    stList *overComps = malnSet_getOverlappingPendingComps(malnSet, joinComp->seq, joinComp->chromStart, joinComp->chromEnd, mafTreeLocRoot, NULL);
    for (int i = 0; i < stList_length(overComps); i++) {
        struct malnComp *overComp = stList_get(overComps, i);
        if ((malnComp_getLoc(overComp) == mafTreeLocRoot) && (overComp->blk != joinBlk) && !overComp->blk->deleted) {
            joinBlk = newJoinBlk = joinBlksAtRoot(malnSet, joinBlk, overComp->blk);
        }
    }
    stList_destruct(overComps);
    return newJoinBlk;
}

/* Join one block with any duplications of that block in the same set. Blocks
 * that are consumed marked as deleted and new blocks added to newBlks */
static void joinBlkWithDups(struct malnSet *malnSet, struct malnBlk *joinBlk, struct malnBlkSet *newBlks) {
    // iterate over each component in joinBlk, looking for overlapping components
    // in other blocks.  A join creates a new block, so we start over until no
    // block is joined with this block.

    bool joinedSome = FALSE;
    struct malnBlk *newJoinBlk;
    do {
        newJoinBlk = joinBlkWithDupsPass(malnSet, joinBlk);
        if (newJoinBlk != NULL) {
            joinBlk = newJoinBlk;
            joinedSome = true;
        }
    } while (newJoinBlk != NULL);
    if (joinedSome) {
        malnBlkSet_add(newBlks, joinBlk);
    }
}

/* Join duplication blocks in a set, which evolver outputs as separate
 * blocks. */
void malnJoinWithinSet_joinDups(struct malnSet *malnSet) {
    struct malnBlkSet *newBlks = malnBlkSet_construct();
    struct malnBlkSetIterator *iter = malnSet_getBlocks(malnSet);
    struct malnBlk *joinBlk;
    while ((joinBlk = malnBlkSetIterator_getNext(iter)) != NULL) {
        if (!joinBlk->deleted) {
            joinBlkWithDups(malnSet, joinBlk, newBlks);
        }
    }
    malnBlkSetIterator_destruct(iter);
    malnSet_deleteDying(malnSet);
    malnSet_addBlks(malnSet, newBlks);
    malnBlkSet_destruct(newBlks);
    malnSet_assert(malnSet);
}

/* compare root components */
static int rootSortCmp(const void *va, const void *vb) {
    struct malnBlk *blkA = (struct malnBlk *)va;
    struct malnBlk *blkB = (struct malnBlk *)vb;
    return malnComp_chromCmp(malnBlk_getRootComp(blkA), malnBlk_getRootComp(blkB));
}

/* build root-component sorted list of blocks */
static stList *buildRootSortedBlks(struct malnSet *malnSet) {
    stList *blks = stList_construct();
    struct malnBlkSetIterator *iter = malnSet_getBlocks(malnSet);
    struct malnBlk *blk;
    while ((blk = malnBlkSetIterator_getNext(iter)) != NULL) {
        stList_append(blks, blk);
    }
    malnBlkSetIterator_destruct(iter);
    stList_sort(blks, rootSortCmp);
    return blks;
}

/* can a component from one block be connected to any other component
 * in the other block in the specified orientation */
static bool canConnectCompBlkOrient(struct malnComp *comp1, struct malnBlk *blk2, int orient) {
    struct malnComp *root2 = malnBlk_getRootComp(blk2);
    for (struct malnComp *comp2 = blk2->comps; comp2 != root2; comp2 = comp2->next) {
        if (malnComp_overlapAdjacentOrient(comp1, comp2, orient)) {
            return true;
        }
    }
    return false;
}

/* can any component from one block be connected to any other component
 * in the other block in the specified orientation */
static bool canConnectBlkBlkOrient(struct malnBlk *blk1, struct malnBlk *blk2, int orient) {
    struct malnComp *root1 = malnBlk_getRootComp(blk1);
    for (struct malnComp *comp1 = blk1->comps; comp1 != root1; comp1 = comp1->next) {
        if (canConnectCompBlkOrient(comp1, blk2, orient)) {
            return true;
        }
    }
    return false;
}

/* can a block be joined with the next one in the list */
static bool canJoinAdjacent(stList *rootSortBlks, int nextIdx) {
    struct malnBlk *blk1 = stList_get(rootSortBlks, nextIdx);
    struct malnComp *root1 = malnBlk_getRootComp(blk1);
    struct malnBlk *blk2 = stList_get(rootSortBlks, nextIdx+1);
    struct malnComp *root2 = malnBlk_getRootComp(blk2);
    if (malnComp_overlap(root1, root2)) {
        return true;  // always join if overlap
    } else if (!malnComp_overlapAdjacent(root1, root2)) {
        return false;  // not overlap or adjacent, can't join roots
    } else {
        int orient = (root1->strand == root2->strand) ? 1 : -1;
        return canConnectBlkBlkOrient(blk1, blk2, orient);
    }
}

/* find a range of blocks to join contiguous in the list, return start of element
 * after the last to join */
static int findAdjacentToJoin(stList *rootSortBlks, int startIdx) {
    int nextIdx = startIdx+1;
    while (nextIdx < stList_length(rootSortBlks) && canJoinAdjacent(rootSortBlks, nextIdx-1)) {
        nextIdx++;
    }
    return nextIdx;
}

/* join a set of adjacent blocks.  return the actually next index, as the number
 * joined might be limited by maxBlkWidth */
static int joinAdjacents(struct malnSet *malnSet, struct malnBlkSet *newBlks, stList *rootSortBlks, int startIdx, int nextIdx, int maxBlkWidth) {
    struct malnBlk *joinedBlk = stList_get(rootSortBlks, startIdx);
    stList_set(rootSortBlks, startIdx, NULL); // just paranoid, don't leave ptrs to blks that might be delete
#if 0 // FIXME:
    malnBlk_dump(joinedBlk, stderr, "joinAdjacents: %d %d", startIdx, nextIdx);
#endif
    int iBlk;
    for (iBlk = startIdx+1; (iBlk < nextIdx) && (joinedBlk->alnWidth <= maxBlkWidth); iBlk++) {
        joinedBlk = joinBlksAtRoot(malnSet, joinedBlk, stList_get(rootSortBlks, iBlk));
        stList_set(rootSortBlks, iBlk, NULL);
    }
    if (joinedBlk->malnSet == NULL) {
        // created this block
        malnBlkSet_add(newBlks, joinedBlk);
    }
    return iBlk;
}

/* join adjacent blocks into a single block, up to the next root gap, or place
 * where there are no components except the root.  Return index of next
 * starting point. */
static int joinSomeAdjacent(struct malnSet *malnSet, struct malnBlkSet *newBlks, stList *rootSortBlks, int startIdx, int maxBlkWidth) {
    int nextIdx = findAdjacentToJoin(rootSortBlks, startIdx);
    nextIdx = joinAdjacents(malnSet, newBlks, rootSortBlks, startIdx, nextIdx, maxBlkWidth);
    return nextIdx;
}

/* Join adjacent and overlapping blocks in a set.  Stop joining at columns
 * where only the root component crosses that column and the root is adjacent
 * and not overlapping. */
void malnJoinWithinSet_joinOverlapAdjacent(struct malnSet *malnSet, int maxBlkWidth) {
    stList *rootSortBlks = buildRootSortedBlks(malnSet);
    struct malnBlkSet *newBlks = malnBlkSet_construct();

    int nextIdx = 0;
    while (nextIdx < stList_length(rootSortBlks)) {
        nextIdx = joinSomeAdjacent(malnSet, newBlks, rootSortBlks, nextIdx, maxBlkWidth);
        malnSet_deleteDying(malnSet);  // keep memory down
    }
    malnSet_addBlks(malnSet, newBlks);
    stList_destruct(rootSortBlks);
    malnBlkSet_destruct(newBlks);
    malnSet_assert(malnSet);
}

