#include "malnJoinSets.h"
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

/* join two components, returning the potentially new joinedBlk */
static struct malnBlk *joinCompComp(struct malnComp *comp1, struct malnComp *comp2, struct malnBlk *joinedBlk) {
    return joinedBlk;
}

/* join for specified joinable component, returning the potentially new joinedBlk */
static struct malnBlk *joinCompSet(struct malnComp *comp1, struct malnSet *malnSet2, struct malnBlk *joinedBlk) {
    stList *overComps = malnSet_getOverlappingPendingComps(malnSet2, comp1->seq, comp1->chromStart, comp1->chromEnd, malnCompTreeRoot|malnCompTreeLeaf);
    // FIXME: struct malnBlk *prevJoinedBlk = joinedBlk;
    for (int i = 0; i < stList_length(overComps); i++) {
        struct malnComp *comp2 = stList_get(overComps, i);
        if (malnComp_canJoin(comp1, comp2)) {
            joinedBlk = joinCompComp(comp1, comp2, joinedBlk);
        }
    }
    stList_destruct(overComps);
    return joinedBlk;
}

/* join blocks overlapping a block of another set. */
static void joinBlkSet(struct malnSet *malnSetJoined, struct malnBlk *blk1, struct malnSet *malnSet2) {
    assert(!blk1->done);
    blk1->done = true;
    struct malnBlk *joinedBlk = NULL;
    for (struct malnComp *comp1 = blk1->comps; comp1 != NULL; comp1 = comp1->next) {
        if (malnComp_joinable(comp1)) {
            joinedBlk = joinCompSet(comp1, malnSet2, joinedBlk);
        }
    }
    if (joinedBlk != NULL) {
        malnSet_addBlk(malnSetJoined, joinedBlk);
    }
}

/* add blocks that were not joined into the alignment */
static void addUndone(struct malnSet *malnSetJoined, struct malnSet *malnSet) {
    stSortedSetIterator *iter = malnSet_getBlocks(malnSet);
    struct malnBlk *blk;
    while ((blk = stSortedSet_getNext(iter)) != NULL) {
        if (!blk->done) {
            malnSet_addBlk(malnSetJoined, malnBlk_constructClone(blk));
            blk->done = true;
        }
    }
    stSortedSet_destructIterator(iter);
}

/* join two sets, generating a third */
struct malnSet *malnJoinSets(struct malnSet *malnSet1, struct malnSet *malnSet2) {
    struct malnSet *malnSetJoined = malnSet_construct(malnSet_getGenomes(malnSet1), malnSet_getRefGenome(malnSet1));

    stSortedSetIterator *iter1 = malnSet_getBlocks(malnSet1);
    struct malnBlk *blk1;
    while ((blk1 = stSortedSet_getNext(iter1)) != NULL) {
        joinBlkSet(malnSetJoined, blk1, malnSet2);
    }
    stSortedSet_destructIterator(iter1);
    addUndone(malnSetJoined, malnSet1);
    addUndone(malnSetJoined, malnSet2);
    return malnSetJoined;
}
