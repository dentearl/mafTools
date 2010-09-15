#include "malnJoinDups.h"
#include "malnSet.h"
#include "malnBlk.h"
#include "malnComp.h"
#include "malnJoinBlks.h"
#include "sonLibSortedSet.h"
#include "sonLibList.h"
#include "sonLibHash.h"
#include "stSafeC.h"
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

/* report a multiple parent, either to stderr or by aborting */
static void reportMultiParent(struct malnComp *joinMultiComp, struct malnComp *dupMultiComp, bool discardTwoParents) {
    char *msg = stSafeCDynFmt("multiple parents detected in components %s:%d-%d (%c) and %s:%d-%d (%c)",
                              joinMultiComp->seq->orgSeqName, joinMultiComp->start, joinMultiComp->end, joinMultiComp->strand,
                              dupMultiComp->seq->orgSeqName, dupMultiComp->start, dupMultiComp->end, dupMultiComp->strand);
    if (discardTwoParents) {
        fprintf(stderr, "Warning: %s\n", msg);
    } else {
        errAbort("Error: %s", msg);
    }
    stSafeCFree(msg);
}

/* Check for multiple parents for non-reference sequence.  Either generate a
 * warning and return true or generate an error. */
static bool checkMultiParent(struct malnComp *joinComp, struct malnComp *dupComp, bool discardTwoParents) {
    for (struct malnComp *comp1 = joinComp->blk->comps; comp1 != NULL; comp1 = comp1->next) {
        if (comp1 != joinComp) {
            for (struct malnComp *comp2 = dupComp->blk->comps; comp2 != NULL; comp2 = comp2->next) {
                if ((comp1 != comp2) && malnComp_overlap(comp1, comp2)) {
                    reportMultiParent(comp1, comp2, discardTwoParents);
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
 * invalidated. */
static void joinCompWithDups(struct malnSet *malnSet, struct malnComp *joinComp, stHash *deleteBlksTbl, bool discardTwoParents) {
    stList *overComps = malnSet_getOverlappingPendingComps(malnSet, joinComp->seq, joinComp->chromStart, joinComp->chromEnd, malnCompTreeRoot);
    struct malnComp *startComp = joinComp;
    int addCnt = 0;
    for (int i = 0; i < stList_length(overComps); i++) {
        struct malnComp *dupComp = stList_get(overComps, i);
        // ignore starting point and deleted
        if ((dupComp != startComp) && (!isDeletedComp(deleteBlksTbl, dupComp))) {
            if (checkMultiParent(joinComp, dupComp, discardTwoParents)) {
                flagDeleted(deleteBlksTbl, dupComp);
            } else {
                joinComp = joinCompWithDup(malnSet, joinComp, dupComp, deleteBlksTbl);
                addCnt++;
            }
        }
    }
    if (addCnt > 0) {
        malnSet_addBlk(malnSet, joinComp->blk);
    }
    stList_destruct(overComps);
}

/* Join duplication blocks in a set, which evolver outputs as separate
 * blocks.  Duplications will only be joined at the root */
void malnJoin_joinSetDups(struct malnSet *malnSet, bool discardTwoParents) {
    stHash *deleteBlksTbl = stHash_construct2(NULL, (void (*)(void *))malnBlk_destruct);
    stSortedSetIterator *iter = malnSet_getBlocks(malnSet);
    struct malnBlk *joinBlk;
    while ((joinBlk = stSortedSet_getNext(iter)) != NULL) {
        if (!isDeletedBlk(deleteBlksTbl, joinBlk)) {
            joinCompWithDups(malnSet, malnBlk_getRootComp(joinBlk), deleteBlksTbl, discardTwoParents);
        }
    }
    stSortedSet_destructIterator(iter);
    stHash_destruct(deleteBlksTbl);
    malnSet_assert(malnSet);
}

