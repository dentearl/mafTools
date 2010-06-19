#include "malnJoinDups.h"
#include "malnSet.h"
#include "malnBlk.h"
#include "malnComp.h"
#include "malnJoinBlks.h"
#include "sonLibSortedSet.h"
#include "sonLibList.h"
#include "genome.h"
#include "common.h"
#include <stdbool.h>
#include <unistd.h>

/* join one into joinComp, returning the new joined comp. */
static struct malnComp *joinCompWithDup(struct malnSet *malnSet, struct malnComp *joinComp, struct malnComp *dupComp, stList *deleteBlkList) {
    assert(!dupComp->blk->done);
    assert(joinComp->blk->done);
    struct malnBlk *joinedBlk = malnJoinBlks(joinComp, dupComp, NULL);
    joinComp->blk->done = true;
    stList_append(deleteBlkList, joinComp->blk);
    dupComp->blk->done = true;
    stList_append(deleteBlkList, dupComp->blk);
    return malnBlk_getRootComp(joinedBlk);
}

/* Join one block with any duplications of that block in the same
 * set.  
 *  - If joinComp doesn't have any dups, it' block is marked as done.
 *  - If joinComp does have dups, it is join with the dups into a
 *    joinedComp. The joinComp and dups are marked as done and added to the
 *    delete list.  The new joinedComp is added to the set and left undone,
 *    so any transitive joins will be picked up on the next round.
 * The delete list is used to delay deletes so that iterators are not
 * invalidated.  Return true if any were joined.
 */
static bool joinCompWithDups(struct malnSet *malnSet, struct malnComp *joinComp, stList *deleteBlkList) {
    joinComp->blk->done = true;  // set upfront so not returned in overlap
    stList *overComps = malnSet_getOverlappingPendingComps(malnSet, joinComp->seq, joinComp->chromStart, joinComp->chromEnd, malnCompTreeRoot);
    for (int i = 0; i < stList_length(overComps); i++) {
        joinComp = joinCompWithDup(malnSet, joinComp, stList_get(overComps, i), deleteBlkList);
    }
    if (overComps != NULL) {
        assert(joinComp->blk->malnSet == NULL);
        malnSet_addBlk(malnSet, joinComp->blk);
    }
    stList_destruct(overComps);
    return (overComps != NULL);
}

/* one pass in duplicate joins, return true if any joined.  Multiple
 * passes are done to handling transitive join cases */
static bool joinSetDupsPass(struct malnSet *malnSet, stList *deleteBlkList) {
    bool joinedSome = false;
    stSortedSetIterator *iter = malnSet_getBlocks(malnSet);
    struct malnBlk *joinBlk;
    while ((joinBlk = stSortedSet_getNext(iter)) != NULL) {
        if (!joinBlk->done) {
            if (joinCompWithDups(malnSet, malnBlk_getRootComp(joinBlk), deleteBlkList)) {
                joinedSome = true;
            }
        }
    }
    stSortedSet_destructIterator(iter);
    return joinedSome;
}

/* Join duplication blocks in a set, which evolver outputs as separate
 * block.  Duplications will only be joined at the root */
void malnJoin_joinSetDups(struct malnSet *malnSet) {
    stList *deleteBlkList = stList_construct3(0, (void (*)(void *))malnBlk_destruct);
    while (joinSetDupsPass(malnSet, deleteBlkList)) {
        continue;
    }
    malnSet_assertDone(malnSet);
    malnSet_clearDone(malnSet);
    stList_destruct(deleteBlkList);
    malnSet_assert(malnSet);
}

