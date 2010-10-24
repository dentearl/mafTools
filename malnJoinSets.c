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
 * Notes:
 * - Initial join happens at the specified join genome, however subsequent joins
 *   can happen on other nodes, due to join genome becoming internal.  In this
 *   case, two root nodes are actually joined into one.  
 */

/* Object use to track the state of blocks that have been process or have been created due
 * to splitting blocks. */
struct blkState {
    struct malnSet *pending1;   // blocks create by splitting blocks from MAF1 that need to be processed
    struct malnSet *pending2;   // blocks create by splitting blocks from MAF2 that need to be processed
    struct malnBlkSet *done2;   // blocks in MAF 2 that have been processed
};
   

/* constructor */
static struct blkState* blkState_construct(struct Genomes *genomes) {
    struct blkState *state;
    AllocVar(state);
    state->pending1 = malnSet_construct(genomes);
    state->done2 = malnBlkSet_construct();
    state->pending2 = malnSet_construct(genomes);
    return state;
}


/* destructor */
static void blkState_destruct(struct blkState* state) {
    malnSet_destruct(state->pending1);
    malnBlkSet_destruct(state->done2);
    malnSet_destruct(state->pending2);
    freeMem(state);
}

/* object used to store current join block/component, which changes when new blocks are
 * merged. */
struct joinBlkComp {
    struct malnBlk *blk;
    struct malnComp *comp;
};

/* split out a region of a block */
static struct malnBlk *splitRegion(struct malnComp *comp, int chromStart, int chromEnd) {
    int start, end, alnStart, alnEnd;
    malnComp_chromRangeToStrandRange(comp, chromStart, chromEnd, &start, &end);
    if (!malnComp_seqRangeToAlnRange(comp, start, end, &alnStart, &alnEnd)) {
        errAbort("BUG: splitRegion barfing");
    }
    return malnBlk_constructSubrange(comp->blk, alnStart, alnEnd);
}

/* split off region of a block and save for future processing */
static void splitOutOfBounds(struct malnComp *belowComp, int chromStart, int chromEnd, struct malnSet *belowPending) {
    struct malnBlk *blk = splitRegion(belowComp, chromStart, chromEnd);
    malnSet_addBlk(belowPending, blk);
}

/* Split out region of a block and return corresponding component. in done set
 * is not NULL, add old block to that set, otherwise delete the old block. Return
 * the equivalent component. */
static struct malnComp *splitInBounds(struct malnComp *belowComp, int chromStart, int chromEnd, struct malnBlkSet *done) {
    assert(malnComp_getLoc(belowComp) == mafTreeLocRoot);
    struct malnBlk *newBlk = splitRegion(belowComp, chromStart, chromEnd);
    struct malnComp *newComp = malnBlk_findCompByChromRange(newBlk, belowComp->seq, chromStart, chromEnd);
    assert(newComp != NULL);

    // dispose of old block
    if (done != NULL) {
        malnBlkSet_add(done, belowComp->blk);
    } else if (belowComp->blk->malnSet != NULL) {
        malnSet_markAsDeleted(belowComp->blk->malnSet, belowComp->blk);
    } else {
        malnBlk_destruct(belowComp->blk);
    }
    return newComp;
}

/* Split a block into up to three parts so that the component being merged
 * from below (component is root) doesn't extend beyond the bounds of the
 * above component (new parent).  The unused ranges of belowBlock are added to
 * the pending set.  If done is not NULL, blocks to no longer use are added to
 * this list, otherwise they are marked for deletion or deleted. No blocks are
 * derived from aboveComp. */
static struct malnComp *splitToBounds(struct malnComp *aboveComp, struct malnComp *belowComp, struct malnSet *belowPending, struct malnBlkSet *done) {
    assert(malnComp_overlap(aboveComp, belowComp));
    assert(malnComp_getLoc(belowComp) == mafTreeLocRoot);
    bool splitSome = false;
    if (belowComp->chromStart < aboveComp->chromStart) {
        // before aboveComp
        splitOutOfBounds(belowComp, belowComp->chromStart, aboveComp->chromStart, belowPending);
        splitSome = true;
    }
    if (belowComp->chromEnd > aboveComp->chromEnd) {
        // after aboveComp
        splitOutOfBounds(belowComp, aboveComp->chromEnd, belowComp->chromEnd, belowPending);
        splitSome = true;
    }
    if (splitSome) {
        return splitInBounds(belowComp, max(belowComp->chromStart, aboveComp->chromStart), min(belowComp->chromEnd, aboveComp->chromEnd), done);
    } else {
        return belowComp;
    }
}

/* join two components updating joining object */
static void joinCompWithComp(struct joinBlkComp *joining, struct malnComp *comp2, struct blkState *state) {
    // trimmed block being attached from below
    if (malnComp_getLoc(joining->comp) == mafTreeLocRoot) {
        joining->comp = splitToBounds(comp2, joining->comp, state->pending1, NULL);
        joining->blk = joining->comp->blk;
    } else {
        comp2 = splitToBounds(joining->comp, comp2, state->pending2, state->done2);
    }

    struct malnComp *joinedComp = NULL;
    struct malnBlk *joinedBlk = malnJoinBlks(joining->comp, comp2, &joinedComp);
    
    // if joining block is not in a set, just free it.
    if (joining->blk->malnSet == NULL) {
        malnBlk_destruct(joining->blk);
    }
    malnBlkSet_add(state->done2, comp2->blk);

    // return resulting block for more additions
    joining->blk = joinedBlk;
    joining->comp = joinedComp;
}

/* join for specified joinable component, returning the potentially new joinedBlk.  Return
 * true if any joined.  This  will be restarted if any join happens  */
static bool joinCompWithSet(struct joinBlkComp *joining, struct malnSet *malnSet2, struct blkState *state) {
    // must stop if we join any, as bounds might be changed by trimming
    bool joinedOne = false;
    stList *overComps2 = malnSet_getOverlappingPendingComps(malnSet2, joining->comp->seq, joining->comp->chromStart, joining->comp->chromEnd, mafTreeLocRoot|mafTreeLocLeaf, state->done2);
    for (int i = 0; i < (stList_length(overComps2) && !joinedOne); i++) {
        struct malnComp *comp2 = stList_get(overComps2, i);
        if ((!malnBlkSet_contains(state->done2, comp2->blk)) && malnComp_canJoin(joining->comp, comp2)) {
            joinCompWithComp(joining, comp2, state);
            joinedOne = true;
        }
    }
    stList_destruct(overComps2);
    return joinedOne;
}

/* one pass over guideGenome components attempt to join with the specified set */
static bool joinBlkWithSetPass( struct joinBlkComp *joining, struct Genome *guideGenome, struct malnBlk *blk1, struct malnSet *malnSet2, struct blkState *state) {
    for (joining->comp = joining->blk->comps; joining->comp != NULL; joining->comp = joining->comp->next) {
        if (joining->comp->seq->genome == guideGenome) {
            if (joinCompWithSet(joining, malnSet2, state)) {
                return true;
            }
        }
    }
    return false;
}

/* join blocks overlapping a block of another set, including any fragments. If none are joined,
 * add blk1, copy blk1 to new set */
static void joinBlkWithSet(struct malnSet *malnSetJoined, struct Genome *guideGenome, struct malnBlk *blk1, struct malnSet *malnSet2, struct blkState *state) {
    struct joinBlkComp joining = {blk1, NULL};
    // Since joining creates a new block and we want to continue to search other components,
    // we start over with the first component when we create a new block.  We do this until
    // we go through all of the joinable components without actually merging a block.
    bool joinedOne = false;
    do {
        joinedOne = joinBlkWithSetPass(&joining, guideGenome, blk1, malnSet2, state);
        if (!joinedOne) {
            joinedOne = joinBlkWithSetPass(&joining, guideGenome, blk1, state->pending2, state);
        }
    } while (joinedOne);
        
    if (joining.blk->malnSet == NULL) {
        // new block, something was joined
        malnSet_addBlk(malnSetJoined, joining.blk);
    } else {
        // nothing joined, copy it.
        assert(joining.blk == blk1);
        malnSet_addBlk(malnSetJoined, malnBlk_constructClone(joining.blk));
    }
}

/* join all pending blocks that have been created by splitting blocks */
static void joinPendingWithSet(struct malnSet *malnSetJoined, struct Genome *guideGenome, struct malnSet *malnSet2, struct blkState *state) {
    // this may create more pending blocks to join
    struct malnBlk *blk1;
    while ((blk1 = malnSet_popBlk(state->pending1)) != NULL) {
        joinBlkWithSet(malnSetJoined, guideGenome, blk1, malnSet2, state);
    }
}

/* join blocks between sets */
static void joinBlksWithSet(struct malnSet *malnSetJoined, struct Genome *guideGenome, struct malnSet *malnSet1, struct malnSet *malnSet2, struct blkState *state) {
    struct malnBlkSetIterator *iter1 = malnSet_getBlocks(malnSet1);
    struct malnBlk *blk1;
    while ((blk1 = malnBlkSetIterator_getNext(iter1)) != NULL) {
        joinBlkWithSet(malnSetJoined, guideGenome, blk1, malnSet2, state);
        joinPendingWithSet(malnSetJoined, guideGenome, malnSet2, state);
    }
    malnBlkSetIterator_destruct(iter1);
}

/* add blocks from second set there  were not joined into the alignment */
static void addUndone(struct malnSet *malnSetJoined, struct malnSet *malnSet2, struct blkState *state) {
    struct malnBlkSetIterator *iter = malnSet_getBlocks(malnSet2);
    struct malnBlk *blk2;
    while ((blk2 = malnBlkSetIterator_getNext(iter)) != NULL) {
        if (!malnBlkSet_contains(state->done2, blk2)) {
            malnSet_addBlk(malnSetJoined, malnBlk_constructClone(blk2));
        }
    }
    malnBlkSetIterator_destruct(iter);
}

/* join two sets, generating a third */
struct malnSet *malnJoinSets(struct Genome *guideGenome, struct malnSet *malnSet1, struct malnSet *malnSet2) {
    struct blkState *state = blkState_construct(malnSet_getGenomes(malnSet1));
    struct malnSet *malnSetJoined = malnSet_construct(malnSet_getGenomes(malnSet1));

    joinBlksWithSet(malnSetJoined, guideGenome, malnSet1, malnSet2, state);
    addUndone(malnSetJoined, malnSet2, state);
    addUndone(malnSetJoined, state->pending2, state);

    assert(malnSet_popBlk(state->pending1) == NULL);
    blkState_destruct(state);
    malnSet_assert(malnSetJoined);
    return malnSetJoined;
}
