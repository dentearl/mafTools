#include "stMalnJoin.h"
#include "stMalnSet.h"
#include "stMalnBlk.h"
#include "stMalnComp.h"
#include "stMalnJoinBlks.h"
#include "genome.h"
#include "common.h"
#include <stdbool.h>
#include <unistd.h>

/* return a reverse list of reference components in a block. */
static struct slRef *getRefComps(struct Genome *refGenome, struct stMalnBlk *blk) {
    struct slRef *refComps = NULL;
    for (struct stMalnComp *comp = blk->comps; comp != NULL; comp = comp->next) {
        if (comp->seq->genome == refGenome) {
            slAddHead(&refComps, slRefNew(comp));
        }
    }
    return refComps;
}

/* recursively find a component overlapping the specified component, or NULL if none found.  Start
 * with the last component and work back to try to get the top of the tree.
 * FIXME: should use the tree. */
static struct stMalnComp *findOverlap(struct stMalnComp *comp1, struct stMalnComp *comps2) {
    struct stMalnComp *over2 = NULL;
    if (comps2 != NULL) {
        over2 = findOverlap(comp1, comps2->next);
        if ((over2 == NULL) && stMalnComp_overlap(comp1, comps2)) {
            over2 = comps2;
        }
    }
    return over2;
}

/* attempt to join one block using a reference component and an overlapping block */
static bool joinBlkCompOverBlk(struct stMalnSet *malnSetJoined, struct stMalnBlk *blk1, struct stMalnComp *refComp1, struct stMalnBlk *blk2) {
    // this should really always find something ...
    struct stMalnComp *refComp2 = findOverlap(refComp1, blk2->comps);
    assert(refComp2 != NULL);
    if (refComp2 != NULL) {
        struct stMalnBlk *blkJoined = stMalnJoinBlks(blk1, refComp1, blk2, refComp2);
        if (blkJoined != NULL) {
            stMalnSet_addBlk(malnSetJoined, blkJoined);
            return true;
        }
    }
    return false;
}

/* attempt to join one block using a reference component */
static bool joinBlkComp(struct stMalnSet *malnSetJoined, struct stMalnBlk *blk1, struct stMalnComp *refComp1, struct stMalnSet *malnSet2) {
    bool got = false;
    struct slRef *overBlks2 = stMalnSet_getOverlapping(malnSet2, refComp1->seq, refComp1->chromStart, refComp1->chromEnd);
    for (struct slRef *overBlkRef2 = overBlks2; (overBlkRef2 != NULL) && (!got); overBlkRef2 = overBlkRef2->next) {
        got = joinBlkCompOverBlk(malnSetJoined, blk1, refComp1, overBlkRef2->val);
    }
    slFreeList(&overBlks2);
    return got;
}

/* join one block */
static void joinBlk(struct Genome *refGenome, struct stMalnSet *malnSetJoined, struct stMalnBlk *blk1, struct stMalnSet *malnSet2) {
    struct slRef *refComps1 = getRefComps(refGenome, blk1);
    for (struct slRef *refComp1Ref = refComps1; refComp1Ref != NULL; refComp1Ref = refComp1Ref->next) {
        if (joinBlkComp(malnSetJoined, blk1, refComp1Ref->val, malnSet2)) {
            break;
        }
    }
    slFreeList(&refComps1);
}

/* add blocks that were not joined into the alignment */
static void addUndone(struct stMalnSet *malnSetJoined, struct stMalnSet *malnSet) {
    for (struct stMalnBlk *blk = stMalnSet_getBlocks(malnSet); blk != NULL; blk = blk->next) {
        if (!blk->done) {
            stMalnSet_addBlk(malnSetJoined, stMalnBlk_constructClone(blk));
            blk->done = true;
        }
    }
}

/* join two sets, generating a third */
struct stMalnSet *stMalnJoin_joinSets(struct Genome *refGenome, struct stMalnSet *malnSet1, struct stMalnSet *malnSet2) {
    struct stMalnSet *malnSetJoined = stMalnSet_construct(stMalnSet_getGenomes(malnSet1), refGenome);
    for (struct stMalnBlk *blk1 = stMalnSet_getBlocks(malnSet1); blk1 != NULL; blk1 = blk1->next) {
        joinBlk(refGenome, malnSetJoined, blk1, malnSet2);
    }
    addUndone(malnSetJoined, malnSet1);
    addUndone(malnSetJoined, malnSet2);
    return malnSetJoined;
}
