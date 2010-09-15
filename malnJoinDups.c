#include "malnJoinDups.h"
#include "malnSet.h"
#include "malnBlk.h"
#include "malnComp.h"
#include "malnJoinBlks.h"
#include "sonLibSortedSet.h"
#include "sonLibList.h"
#include "sonLibHash.h"
#include "genome.h"
#include "common.h"
#include <stdbool.h>
#include <unistd.h>

/* add a component's block to the delete table */
static void flagDeleted(stHash *deleteBlksTbl, struct malnComp *comp) {
    stHash_insert(deleteBlksTbl, comp->blk, comp->blk);
}

/* is a block in the delete table? */
static bool isDeletedBlk(stHash *deleteBlksTbl, struct malnBlk *blk) {
    return stHash_search(deleteBlksTbl, blk) != NULL; 
}

/* is a component in the delete table? */
static bool isDeletedComp(stHash *deleteBlksTbl, struct malnComp *comp) {
    return isDeletedBlk(deleteBlksTbl, comp->blk);
}


/* FIXME hack: malnJoinBlks was written assuming none of the non-reference
 * blocks overlap.  Then the code was misused to join tandem dups in the
 * same alignment set.  This broke things badly and means some major work
 * in malnJoinBlks, so we introduced this hack.
 * FIXME lots of debugging code just for this in malnJoinBlks.c. */
static bool isNastyTandemDupCase(struct malnComp *joinComp, struct malnComp *dupComp) {
    if (false) { // set for debugging
        return false; 
    }
    for (struct malnComp *comp1 = joinComp->blk->comps; comp1 != NULL; comp1 = comp1->next) {
        if (comp1 != joinComp) {
            for (struct malnComp *comp2 = dupComp->blk->comps; comp2 != NULL; comp2 = comp2->next) {
                if ((comp1 != comp2) && malnComp_overlap(comp1, comp2)) {
                    return true;
                }
            }
        }
    }
    return false;
}

/* join two blocks associated components, create and third block. Return that
 * blocks root component. Joined blocks are inserted into delete table */
static struct malnComp *joinCompWithDup(struct malnSet *malnSet, struct malnComp *comp1, struct malnComp *comp2, stHash *deleteBlksTbl) {
    struct malnBlk *joinedBlk = malnJoinBlks(comp1, comp2, NULL);
    flagDeleted(deleteBlksTbl, comp1);
    flagDeleted(deleteBlksTbl, comp2);
    return malnBlk_getRootComp(joinedBlk);
}

/* Join one block with any duplications of that block in the same
 * set. Duplicates are added to a table and skipped so that iterators are not
 * invalidated.
 */
static void joinCompWithDups(struct malnSet *malnSet, struct malnComp *joinComp, stHash *deleteBlksTbl) {
    stList *overComps = malnSet_getOverlappingPendingComps(malnSet, joinComp->seq, joinComp->chromStart, joinComp->chromEnd, malnCompTreeRoot);
    struct malnComp *startComp = joinComp;
    int addCnt = 0;
    for (int i = 0; i < stList_length(overComps); i++) {
        struct malnComp *dupComp = stList_get(overComps, i);
        // ignore starting point and deleted
        if ((dupComp != startComp) && (!isDeletedComp(deleteBlksTbl, dupComp)) && !isNastyTandemDupCase(joinComp, dupComp)) {
            joinComp = joinCompWithDup(malnSet, joinComp, dupComp, deleteBlksTbl);
            addCnt++;
        }
    }
    if (addCnt > 0) {
        malnSet_addBlk(malnSet, joinComp->blk);
    }
    stList_destruct(overComps);
}

/* Join duplication blocks in a set, which evolver outputs as separate
 * blocks.  Duplications will only be joined at the root */
void malnJoin_joinSetDups(struct malnSet *malnSet) {
    stHash *deleteBlksTbl = stHash_construct2(NULL, (void (*)(void *))malnBlk_destruct);
    stSortedSetIterator *iter = malnSet_getBlocks(malnSet);
    struct malnBlk *joinBlk;
    while ((joinBlk = stSortedSet_getNext(iter)) != NULL) {
        if (!isDeletedBlk(deleteBlksTbl, joinBlk)) {
            joinCompWithDups(malnSet, malnBlk_getRootComp(joinBlk), deleteBlksTbl);
        }
    }
    stSortedSet_destructIterator(iter);
    stHash_destruct(deleteBlksTbl);
    malnSet_assert(malnSet);
}

