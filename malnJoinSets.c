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

static bool debug = false;  // FIXME: tmp

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
static void splitOutOfBounds(struct malnComp *belowComp, int chromStart, int chromEnd, struct malnBlkSet *pending) {
    struct malnBlk *blk = splitRegion(belowComp, chromStart, chromEnd);
    malnBlkSet_add(pending, blk);
    if (debug) {
        malnBlk_dump(blk, stderr, "add to pending");
    }
}

/* Split out region of a block and return corresponding component. If in done
 * is not NULL, add old block to that set, otherwise delete the old
 * block. Return the equivalent component. */
static struct malnComp *splitInBounds(struct malnComp *belowComp, int chromStart, int chromEnd) {
    assert(malnComp_getLoc(belowComp) == mafTreeLocRoot);
    struct malnBlk *newBlk = splitRegion(belowComp, chromStart, chromEnd);
    struct malnComp *newComp = malnBlk_findCompByChromRange(newBlk, belowComp->seq, chromStart, chromEnd);
    assert(newComp != NULL);
    malnBlk_markOrDelete(belowComp->blk);
    return newComp;
}

/* Split a block into up to three parts so that the component being merged
 * from below (component is root) doesn't extend beyond the bounds of the
 * above component (new parent).  The unused ranges of belowBlock are added to
 * the pending set.  If done is not NULL, blocks to no longer use are added to
 * this list, otherwise they are marked for deletion or deleted. No blocks are
 * derived from aboveComp. */
static struct malnComp *splitToBounds(struct malnComp *aboveComp, struct malnComp *belowComp, struct malnBlkSet *pending) {
    assert(malnComp_overlap(aboveComp, belowComp));
    assert(malnComp_getLoc(belowComp) == mafTreeLocRoot);
    bool splitSome = false;
    if (belowComp->chromStart < aboveComp->chromStart) {
        // before aboveComp
        splitOutOfBounds(belowComp, belowComp->chromStart, aboveComp->chromStart, pending);
        splitSome = true;
    }
    if (belowComp->chromEnd > aboveComp->chromEnd) {
        // after aboveComp
        splitOutOfBounds(belowComp, aboveComp->chromEnd, belowComp->chromEnd, pending);
        splitSome = true;
    }
    if (splitSome) {
        return splitInBounds(belowComp, max(belowComp->chromStart, aboveComp->chromStart), min(belowComp->chromEnd, aboveComp->chromEnd));
    } else {
        return belowComp;
    }
}

/* join two blocks at the specified components, returning resulting block */
static struct malnBlk *joinCompWithComp(struct malnComp *comp1, struct malnComp *comp2, struct malnBlkSet *pending) {
    if (debug) {
        malnBlk_dump(comp1->blk, stderr, "blk1 before split");
        malnBlk_dump(comp2->blk, stderr, "blk2 before split");
    }
    assert(comp1 != comp2);
    // trimmed block being attached from below
    if (malnComp_getLoc(comp1) == mafTreeLocRoot) {
        comp1 = splitToBounds(comp2, comp1, pending);
    } else {
        comp2 = splitToBounds(comp1, comp2, pending);
    }

    if (debug) {
        malnBlk_dump(comp1->blk, stderr, "blk1 after split");
        malnBlk_dump(comp2->blk, stderr, "blk2 after split");
    }
    struct malnBlk *newBlk = malnJoinBlks(comp1, comp2);
    malnBlk_markOrDelete(comp1->blk);
    malnBlk_markOrDelete(comp2->blk);
    return newBlk;
}

/* join for specified joinable component, returning the potentially new joinedBlk.  Return
 * new joinBlk if any joined, otherwise NULL.  */
static struct malnBlk *joinCompWithSet(struct malnComp *joinComp, struct malnSet *targetSet, struct malnBlkSet *pending) {
    if (debug) {
        malnComp_dump(joinComp, stderr, "joinCompWithSet");
        malnBlk_dump(joinComp->blk, stderr, "  joinComp->blk");
    }
    // must stop if we join any, as bounds might be changed by trimming
    struct malnBlk *newJoinBlk = NULL;
    stList *overComps = malnSet_getOverlappingPendingComps(targetSet, joinComp->seq, joinComp->chromStart, joinComp->chromEnd, mafTreeLocAll, NULL);
    for (int i = 0; i < (stList_length(overComps) && (newJoinBlk == NULL)); i++) {
        struct malnComp *targetComp = stList_get(overComps, i);
        if (debug) {
            malnComp_dump(targetComp, stderr, " overlap:");
        }
        if (malnComp_canJoin(joinComp, targetComp)) {
            if (debug) {
                malnBlk_dump(targetComp->blk, stderr, " joinwith:");
            }
            newJoinBlk = joinCompWithComp(joinComp, targetComp, pending);
            if (debug) {
                malnBlk_dump(newJoinBlk, stderr, " newBlk:");
            }
        }
    }
    stList_destruct(overComps);
    return newJoinBlk;
}

/* One pass over guideGenome components attempt to join with the specified set.  Return new joinBlk
 * if joined, otherwise return NULL if not changed. */
static struct malnBlk *joinBlkWithSetPass(struct Genome *guideGenome, struct malnBlk *joinBlk, struct malnSet *targetSet, struct malnBlkSet *pending) {
    for (struct malnComp *joinComp = joinBlk->comps; joinComp != NULL; joinComp = joinComp->next) {
        if (joinComp->seq->genome == guideGenome) {
            struct malnBlk *newJoinBlk = joinCompWithSet(joinComp, targetSet, pending);
            if (newJoinBlk != NULL) {
                return newJoinBlk;
            }
        }
    }
    return NULL;
}

/* add pending blocks back to an arbitrary set */
static void addPending(struct malnBlkSet *pending, struct malnSet *malnSet) {
    struct malnBlk *blk;
    while ((blk = malnBlkSet_pop(pending)) != NULL) {
        malnSet_addBlk(malnSet, blk);
    }
}

/* join a join a blk with overlapping components, making multiple passes until no more blocks can be
 * joined with it */
static bool joinBlkWithSet(struct malnSet *malnSetJoined, struct Genome *guideGenome, struct malnBlk *blk, struct malnSet *malnSet1, struct malnSet *malnSet2, struct malnBlkSet *pending) {
    bool anyJoined = false;
    // make copy and delete from set so it's not overlapped
    struct malnBlk *joinBlk = malnBlk_constructClone(blk);
    malnBlk_markOrDelete(blk);
    struct malnBlk *newJoinBlk = NULL;
    do {
        // try both ways due to transitive joins across the sets
        newJoinBlk = joinBlkWithSetPass(guideGenome, joinBlk, malnSet2, pending);
        if (newJoinBlk == NULL) {
            newJoinBlk = joinBlkWithSetPass(guideGenome, joinBlk, malnSet1, pending);
        }
        if (newJoinBlk != NULL) {
            addPending(pending, malnSet1);
            joinBlk = newJoinBlk;
            anyJoined = true;
        }
    } while (newJoinBlk != NULL);
    malnSet_addBlk(malnSetJoined, joinBlk);
    return anyJoined;
}

/* join blocks from one set */
static bool joinBlksWithSet(struct malnSet *malnSetJoined, struct Genome *guideGenome, struct malnSet *malnSet1, struct malnSet *malnSet2) {
    // make copy of container for blocks so underlying set can be modified
    // during loop (doesn't copy actual blocks)
    bool anyJoined = false;
    struct malnBlkSet *pending = malnBlkSet_construct();
    struct malnBlkSet *set1Blks = malnSet_getBlockSetCopy(malnSet1);
    struct malnBlkSetIterator *iter1 = malnBlkSet_getIterator(set1Blks);
    struct malnBlk *nextBlk;
    while ((nextBlk = malnBlkSetIterator_getNext(iter1)) != NULL) {
        // make copy and delete from set so it's not overlapped
        if (joinBlkWithSet(malnSetJoined, guideGenome, nextBlk, malnSet1, malnSet2, pending)) {
            anyJoined = true;
        }
    }
    malnBlkSetIterator_destruct(iter1);
    malnBlkSet_destruct(set1Blks);
    malnBlkSet_destruct(pending);
    return anyJoined;
}

/* join blocks from both sets with both sets, continuing until nothing done  */
static void joinBlksWithSets(struct malnSet *malnSetJoined, struct Genome *guideGenome, struct malnSet *malnSet1, struct malnSet *malnSet2) {
    bool anyJoined;
    do {
        anyJoined = false;
        if (joinBlksWithSet(malnSetJoined, guideGenome, malnSet1, malnSet2)) {
            anyJoined = true;
        }
        if (joinBlksWithSet(malnSetJoined, guideGenome, malnSet2, malnSet1)) {
            anyJoined = true;
        }
    } while (anyJoined);
}

/* add blocks from set that were not joined into the alignment */
static void addUndone(struct malnSet *malnSetJoined, struct malnSet *malnSet) {
    // FIXME: still needed??
    // ones marked deleted are automatically skipped
    struct malnBlkSetIterator *iter = malnSet_getBlocks(malnSet);
    struct malnBlk *blk;
    while ((blk = malnBlkSetIterator_getNext(iter)) != NULL) {
        malnSet_addBlk(malnSetJoined, malnBlk_constructClone(blk));
    }
    malnBlkSetIterator_destruct(iter);
}

/* join two sets, generating a third */
struct malnSet *malnJoinSets(struct Genome *guideGenome, struct malnSet *malnSet1, struct malnSet *malnSet2) {
    struct malnSet *malnSetJoined = malnSet_construct(malnSet_getGenomes(malnSet1), NULL);
    joinBlksWithSets(malnSetJoined, guideGenome, malnSet1, malnSet2);
    addUndone(malnSetJoined, malnSet1);
    addUndone(malnSetJoined, malnSet2);
    malnSet_assert(malnSetJoined);
    return malnSetJoined;
}
