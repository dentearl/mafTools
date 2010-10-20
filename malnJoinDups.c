#include "malnJoinDups.h"
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

static const bool debug = false; // FIXME: tmp

/* join two blocks by associated root components, create and third
 * block. Return that new blocks root component. Joined blocks are inserted into
 * delete table */
static struct malnBlk *joinCompWithDup(struct malnSet *malnSet, struct malnBlk *blk1, struct malnBlk *blk2) {
    struct malnBlk *joinBlk = malnJoinBlks(malnBlk_getRootComp(blk1), malnBlk_getRootComp(blk2), NULL);
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
        if (malnComp_canJoin(joinComp, overComp) && (overComp->blk != joinBlk) && !overComp->blk->deleted) {
            joinBlk = newJoinBlk = joinCompWithDup(malnSet, joinBlk, overComp->blk);
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
void malnJoinDups_joinSetDups(struct malnSet *malnSet) {
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

