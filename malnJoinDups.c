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

/* join two blocks associated components, create and third block. Return that
 * blocks root component. Joined blocks are inserted into delete table */
static struct malnComp *joinCompWithDup(struct malnSet *malnSet, struct malnComp *comp1, struct malnComp *comp2) {
    if (debug) {
        malnComp_dump(comp2, "joinCompWithDup comp2", stderr);
    }
    struct malnBlk *joinedBlk = malnJoinBlks(comp1, comp2, NULL);
    malnBlk_markOrDelete(comp1->blk);
    malnBlk_markOrDelete(comp2->blk);
    return malnBlk_getRootComp(joinedBlk);
}

/* attempt to join a block with other blocks using the specified component.
 * Return updated block when one join is achieved, or NULL if no join was done. */
static struct malnBlk *joinCompWithDups(struct malnSet *malnSet, struct malnBlk *joinBlk, struct malnComp *joinComp, struct malnBlkSet *doneBlks) {
    if (debug) {
        malnComp_dump(joinComp, "joinCompWithDups", stderr);
    }
    malnBlkSet_add(doneBlks, joinComp->blk);

    stList *overComps = malnSet_getOverlappingAdjacentPendingComps(malnSet, joinComp->seq, joinComp->chromStart, joinComp->chromEnd, mafTreeLocAll, doneBlks);
    for (int i = 0; i < stList_length(overComps); i++) {
        struct malnComp *dupComp = stList_get(overComps, i);
        if (malnComp_canJoin(joinComp, dupComp) && !dupComp->blk->deleted) {
            joinComp = joinCompWithDup(malnSet, joinComp, dupComp);
            malnBlkSet_add(doneBlks, joinComp->blk);
        }
    }
    stList_destruct(overComps);
    return (joinComp->blk != joinBlk) ? joinComp->blk : NULL;
}

/* Join one block with any duplications of that block in the same
 * set. Duplicates are added to a table and skipped so that iterators are not
 * invalidated. */
static void joinBlkWithDups(struct malnSet *malnSet, struct malnBlk *joinBlk, struct malnBlkSet *doneBlks) {
    // iterate over each component in joinBlk, looking for overlapping components
    // in other blocks.  A join creates a new block, so we start over until no
    // block is joined with this block.

    if (debug) { // FIXME: tmp
        malnBlk_dump(joinBlk, "joinBlkWithDups", stderr);
    }
    bool joinedSome = FALSE;
    struct malnBlk *newJoinBlk;
    do {
        newJoinBlk = NULL;
        for (struct malnComp *joinComp = joinBlk->comps; (joinComp != NULL) && (newJoinBlk == NULL); joinComp = joinComp->next) {
            newJoinBlk = joinCompWithDups(malnSet, joinBlk, joinComp, doneBlks);
        }
        if (newJoinBlk != NULL) {
            joinBlk = newJoinBlk;
            joinedSome = true;
        }
    } while (newJoinBlk != NULL);
    if (joinedSome) {
        malnSet_addBlk(malnSet, joinBlk);
    }
}

/* Join duplication blocks in a set, which evolver outputs as separate
 * blocks. */
void malnJoinDups_joinSetDups(struct malnSet *malnSet) {
    struct malnBlkSet *doneBlks = malnBlkSet_construct();
    struct malnBlkSetIterator *iter = malnSet_getBlocks(malnSet);
    struct malnBlk *joinBlk;
    while ((joinBlk = malnBlkSetIterator_getNext(iter)) != NULL) {
        if (!joinBlk->deleted) {
            joinBlkWithDups(malnSet, joinBlk, doneBlks);
        }
    }
    malnBlkSetIterator_destruct(iter);
    malnSet_deleteDying(malnSet);
    malnSet_assert(malnSet);
    malnBlkSet_destruct(doneBlks);
}

