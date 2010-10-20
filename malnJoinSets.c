#include "malnJoinSets.h"
#include "malnSet.h"
#include "malnBlk.h"
#include "malnBlkSet.h"
#include "malnComp.h"
#include "malnJoinBlks.h"
#include "sonLibList.h"
#include "genome.h"
#include "common.h"
#include <stdbool.h>
#include <unistd.h>

/*
 * Initial join happens at the specified join genome, however subsequent joins
 * can happen on other nodes, due to join genome becoming internal.  In this
 * case, two root nodes are actually joined into one.  
 * This is a* hack around the join not being tree based and the merge not happening
 * as part of the join.
 * 
 *
 * join on sHuman-sChimp
 * initial:
 *   set a:
 *     a1  sHuman-sChimp.chr20 743840-754280 [10440] root
 *         simHuman.chr20      757705-757749    [44] leaf
 * 
 *   set b
 *     b1  sG-sH-sC.chr20      741388-751692 [10304] root
 *         sHuman-sChimp.chr20 743840-754142 [10302] leaf
 * 
 *     b2  sG-sH-sC.chr20      751693-751830 [137] root
 *         sHuman-sChimp.chr20 754143-754280 [137] leaf
 * 
 * 
 * merge1: 
 *    ab1 sG-sH-sC.chr20      741388-751692 [10304] root
 *        sHuman-sChimp.chr20 743840-754280 [10440] internal     
 *        simHuman.chr20      757705-757749    [44] leaf
 * 
 *     b2 sG-sH-sC.chr20      751693-751830 [137] root
 *        sHuman-sChimp.chr20 754143-754280 [137] leaf
 * 
 * merge2: 
 *    ab1 sG-sH-sC.chr20      741388-751692 [10304] root
 *        sHuman-sChimp.chr20 743840-754280 [10440] internal     
 *        simHuman.chr20      757705-757749    [44] leaf
 * 
 *     b2 sG-sH-sC.chr20      751693-751830 [137] root
 *        sHuman-sChimp.chr20 754143-754280 [137] leaf
 */

static bool debug = false;  // FIXME: tmp

/* object used to store current join block/component, which changes when new blocks are
 * merged. */
struct joinBlkComp {
    struct malnBlk *blk;
    struct malnComp *comp;
};

/* join two components updating joining object */
static void joinCompWithComp(struct joinBlkComp *joining, struct malnComp *comp2, struct malnBlkSet *doneBlks) {
    struct malnComp *joinedComp = NULL;
    struct malnBlk *joinedBlk = malnJoinBlks(joining->comp, comp2, &joinedComp);
    
    // if joining block is an a set, then just mark as done, otherwise free it.
    if (joining->blk->malnSet != NULL) {
        malnBlkSet_add(doneBlks, joining->blk);
    } else {
        malnBlk_destruct(joining->blk);
    }
    malnBlkSet_add(doneBlks, comp2->blk);
    // return resulting block for more additions
    joining->blk = joinedBlk;
    joining->comp = joinedComp;
}

/* join for specified joinable component, returning the potentially new joinedBlk.  Return
 * true if any joined. */
static bool joinCompWithSet(struct joinBlkComp *joining, struct malnSet *malnSet2, struct malnBlkSet *doneBlks) {
    bool joinedOne = false;
    stList *overComps2 = malnSet_getOverlappingAdjacentPendingComps(malnSet2, joining->comp->seq, joining->comp->chromStart, joining->comp->chromEnd, mafTreeLocRoot|mafTreeLocLeaf, doneBlks);
    for (int i = 0; i < stList_length(overComps2); i++) {
        struct malnComp *comp2 = stList_get(overComps2, i);
        if (debug) {
            malnComp_dump(joining->comp, stderr, "OVER");
            fprintf(stderr, "\t%d & %d: ", (!malnBlkSet_contains(doneBlks, comp2->blk)), malnComp_canJoin(joining->comp, comp2)); malnComp_dump(comp2, stderr, "comp2");
        }
        if ((!malnBlkSet_contains(doneBlks, comp2->blk)) && malnComp_canJoin(joining->comp, comp2)) {
            joinCompWithComp(joining, comp2, doneBlks);
            joinedOne = true;
        }
    }
    if (debug) {
        fprintf(stderr, "joinedOne: %s: =========================\n", (joinedOne ?"YES":"NO"));
    }
    stList_destruct(overComps2);
    return joinedOne;
}

/* join blocks overlapping a block of another set. If none are joined,
 * add blk1, copy blk1 to new set */
static void joinBlkWithSet(struct malnSet *malnSetJoined, struct Genome *refGenome, struct malnBlk *blk1, struct malnSet *malnSet2, int maxBlkWidth, struct malnBlkSet *doneBlks) {
    struct joinBlkComp joining = {blk1, NULL};
    // since joining creates a new block and we want to continue to search other components,
    // we start over with the first component when we create a new block.  We do this until
    // we go through all of the joinable components without actually merging a block.
    bool joinedOne = false;
    do {
        joinedOne = false;
        for (joining.comp = joining.blk->comps; joining.comp != NULL; joining.comp = joining.comp->next) {
            if (joining.comp->seq->genome == refGenome) {
                if (joinCompWithSet(&joining, malnSet2, doneBlks)) {
                    joinedOne = true;
                }
            }
        }
    } while (joinedOne && (joining.blk->alnWidth < maxBlkWidth));
        
    if (joining.blk->malnSet == NULL) {
        // new block, something was joined
        malnSet_addBlk(malnSetJoined, joining.blk);
    } else {
        // nothing joined, copy it.
        assert(joining.blk == blk1);
        malnSet_addBlk(malnSetJoined, malnBlk_constructClone(joining.blk));
    }
}

/* add blocks that were not joined into the alignment */
static void addUndone(struct malnSet *malnSetJoined, struct malnSet *malnSet, struct malnBlkSet *doneBlks) {
    struct malnBlkSetIterator *iter = malnSet_getBlocks(malnSet);
    struct malnBlk *blk;
    while ((blk = malnBlkSetIterator_getNext(iter)) != NULL) {
        if (!malnBlkSet_contains(doneBlks, blk)) {
            malnSet_addBlk(malnSetJoined, malnBlk_constructClone(blk));
            malnBlkSet_add(doneBlks, blk);
        }
    }
    malnBlkSetIterator_destruct(iter);
}

/* join two sets, generating a third */
struct malnSet *malnJoinSets(struct Genome *refGenome, struct malnSet *malnSet1, struct malnSet *malnSet2, int maxBlkWidth) {
    struct malnBlkSet *doneBlks = malnBlkSet_construct();
    struct malnSet *malnSetJoined = malnSet_construct(malnSet_getGenomes(malnSet1));
    struct malnBlkSetIterator *iter1 = malnSet_getBlocks(malnSet1);
    struct malnBlk *blk1;
    while ((blk1 = malnBlkSetIterator_getNext(iter1)) != NULL) {
        if (!malnBlkSet_contains(doneBlks, blk1)) {
            joinBlkWithSet(malnSetJoined, refGenome, blk1, malnSet2, maxBlkWidth, doneBlks);
        }
    }
    malnBlkSetIterator_destruct(iter1);
    addUndone(malnSetJoined, malnSet2, doneBlks);
    malnBlkSet_destruct(doneBlks);
    malnSet_assert(malnSetJoined);
    return malnSetJoined;
}
