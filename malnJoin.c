#include "malnJoin.h"
#include "malnSet.h"
#include "malnBlk.h"
#include "malnComp.h"
#include "malnJoinBlks.h"
#include "sonLibSortedSet.h"
#include "genome.h"
#include "common.h"
#include <stdbool.h>
#include <unistd.h>

/* join duplication blocks in a set, which */
void malnJoin_joinSetDups(struct malnSet *malnSet) {
    
}

#if 0
/* get the component that can be joined with the specified component, or null
 * if none */
static struct malnComp *getJoinComp(struct malnComp *refComp, struct malnBlk *blk, unsigned treeLocFilter) {
    for (struct malnComp *comp = blk->comps; comp != NULL; comp = comp->next) {
        if (comp->isReference && ((comp->treeLoc & treeLocFilter) != 0) && malnComp_overlap(comp, refComp)) {
            return comp;
        }
    }
    return NULL;
}


/* Add blocks overlapping the particular reference component to the overlapping set */
static void findJoinCompOverlapBlks(struct malnComp *refComp, struct malnSet *malnSet, unsigned treeLocFilter, stSortedSet *overlapping) {
    struct slRef *overBlks = malnSet_getOverlapping(malnSet, refComp->seq, refComp->chromStart, refComp->chromEnd);
    for (struct slRef *overBlkRef = overBlks; overBlkRef != NULL; overBlkRef = overBlkRef->next) {
        struct malnBlk *overBlk = overBlkRef->val;
        if ((!overBlk->done) && (getJoinComp(refComp, overBlk, treeLocFilter) != NULL)) {
            stSortedSet_insert(overlapping, overBlk);
        }
    }
    slFreeList(&overBlks2);
}

/* get a set of blocks with join-reference components overlapping the reference in the
 * seed block. */
static stSortedSet *findJoinOverlap(struct malnBlk *refBlk, struct malnSet *malnSet, unsigned treeLocFilter) {
    stSortedSet *overlapping = stSortedSet_construct();
    for (struct malnComp *refComp = refBlk->comps; refComp != NULL; refComp = refComp->next) {
        findJoinCompOverlapBlks(refComp, malnSet, treeLocFilter, overlapping);
    }
    return overlapping;
}    

/* recursively find a component overlapping the specified component, or NULL if none found.  Start
 * with the last component and work back to try to get the top of the tree.
 * FIXME: should use the tree. */
static struct malnComp *findOverlap(struct malnComp *comp1, struct malnComp *comps2) {
    struct malnComp *over2 = NULL;
    if (comps2 != NULL) {
        over2 = findOverlap(comp1, comps2->next);
        if ((over2 == NULL) && malnComp_overlap(comp1, comps2)) {
            over2 = comps2;
        }
    }
    return over2;
}

/* attempt to join one block using a reference component and an overlapping block */
static bool joinBlkCompOverBlk(struct malnSet *malnSetJoined, struct malnBlk *blk1, struct malnComp *refComp1, struct malnBlk *blk2) {
    // this should really always find something ...
    struct malnComp *refComp2 = findOverlap(refComp1, blk2->comps);
    assert(refComp2 != NULL);
    if (refComp2 != NULL) {
        struct malnBlk *blkJoined = malnJoinBlks(blk1, refComp1, blk2, refComp2);
        if (blkJoined != NULL) {
            malnSet_addBlk(malnSetJoined, blkJoined);
            return true;
        }
    }
    return false;
}

/* attempt to join one block using a reference component */
static bool joinBlkComp(struct malnSet *malnSetJoined, struct malnBlk *blk1, struct malnComp *refComp1, struct malnSet *malnSet2) {
    bool got = false;
    struct slRef *overBlks2 = malnSet_getOverlapping(malnSet2, refComp1->seq, refComp1->chromStart, refComp1->chromEnd);
    for (struct slRef *overBlkRef2 = overBlks2; (overBlkRef2 != NULL) && (!got); overBlkRef2 = overBlkRef2->next) {
        got = joinBlkCompOverBlk(malnSetJoined, blk1, refComp1, overBlkRef2->val);
    }
    slFreeList(&overBlks2);
    return got;
}

/* join one block */
static void joinBlk(struct Genome *refGenome, struct malnSet *malnSetJoined, struct malnBlk *blk1, struct malnSet *malnSet2) {
    struct slRef *refComps1 = getRefComps(refGenome, blk1);
    for (struct slRef *refComp1Ref = refComps1; refComp1Ref != NULL; refComp1Ref = refComp1Ref->next) {
        if (joinBlkComp(malnSetJoined, blk1, refComp1Ref->val, malnSet2)) {
            break;
        }
    }
    slFreeList(&refComps1);
}

/* join on one seed block range from both sets */
static void joinBlk(struct Genome *refGenome, struct malnSet *malnSetJoined, struct malnBlk *blk1, struct malnSet *malnSet2) {
    struct slRef *refComps1 = getRefComps(malnSet_getRefGenome(malnSetJoined), blk1);
    for (struct slRef *refComp1Ref = refComps1; refComp1Ref != NULL; refComp1Ref = refComp1Ref->next) {
        if (joinBlkComp(malnSetJoined, blk1, refComp1Ref->val, malnSet2)) {
            break;
        }
    }
    slFreeList(&refComps1);
}

/* join from two sets based one reference components in a seed block range */
static void joinFromSeedBlk(struct malnSet *malnSetJoined, struct malnBlk *seedBlk, struct malnSet *malnSet1, struct malnSet *malnSet2) {
    static stSortedSet *joinBlks1 = findJoinOverlap(seedBlk, malnSet1);
    static stSortedSet *joinBlks2 = findJoinOverlap(seedBlk, malnSet2);
    

    struct slRef *refComps1 = getRefComps(malnSet_getRefGenome(malnSetJoined), blk1);
    for (struct slRef *refComp1Ref = refComps1; refComp1Ref != NULL; refComp1Ref = refComp1Ref->next) {
        if (joinBlkComp(malnSetJoined, blk1, refComp1Ref->val, malnSet2)) {
            break;
        }
    }
    slFreeList(&refComps1);
    stSortedSet_destruct(joinBlks1);
    stSortedSet_destruct(joinBlks2);
}

/* add blocks that were not joined into the alignment */
static void addUndone(struct malnSet *malnSetJoined, struct malnSet *malnSet) {
    for (struct malnBlk *blk = malnSet_getBlocks(malnSet); blk != NULL; blk = blk->next) {
        if (!blk->done) {
            assert(false); // FIXME: don't this should happen any more
            malnSet_addBlk(malnSetJoined, malnBlk_constructClone(blk));
            blk->done = true;
        }
    }
}

/* join two sets, generating a third */
struct malnSet *malnJoin_joinSets(struct malnSet *malnSet1, struct malnSet *malnSet2) {
    struct malnSet *malnSetJoined = malnSet_construct(malnSet_getGenomes(malnSet1), malnSet_getRefGenome(malnSet1));
    for (struct malnBlk *seedBlk = malnSet_getBlocks(malnSet); seedBlk != NULL; seedBlk1 = seedBlk->next) {
        if (!seedBlk->done) {
            joinFromSeedBlk(refGenome, malnSetJoined, seedblk, malnSet1, malnSet2);
        }
    }
    addUndone(malnSetJoined, malnSet1);
    addUndone(malnSetJoined, malnSet2);
    return malnSetJoined;
}
#endif
