#include "malnMultiParents.h"
#include "malnSet.h"
#include "malnBlk.h"
#include "malnComp.h"
#include "malnDeleteBlks.h"
#include "sonLibSortedSet.h"
#include "sonLibList.h"
#include "stSafeC.h"

/* report a multiple parent, either to stderr or by aborting */
static void reportMultiParent(struct malnComp *comp, struct malnComp *overComp, bool discardTwoParents) {
    char *msg = stSafeCDynFmt("multiple parents detected in components %s:%d-%d (%c) and %s:%d-%d (%c)",
                              comp->seq->orgSeqName, comp->start, comp->end, comp->strand,
                              overComp->seq->orgSeqName, overComp->start, overComp->end, overComp->strand);
    if (discardTwoParents) {
        fprintf(stderr, "Warning: %s\n", msg);
    } else {
        errAbort("Error: %s", msg);
    }
    stSafeCFree(msg);
}

/* Check for multiple parents for non-reference sequence.  Either generate a
 * warning and return true or generate an error. */
static bool checkMultiParent(struct malnComp *comp, struct malnComp *overComp, bool discardTwoParents) {
    for (struct malnComp *comp2 = overComp->blk->comps; comp2 != NULL; comp2 = comp2->next) {
        if ((comp != comp2) && malnComp_overlap(comp, comp2)) {
            reportMultiParent(comp, comp2, discardTwoParents);
            return true;
        }
    }
    return false;
}

/* compare a component to other components */
static void checkCompForMultiParents(struct malnSet *malnSet, struct malnComp *comp, struct malnDeleteBlks *delBlks, bool discardTwoParents) {
    stList *overComps = malnSet_getOverlappingPendingComps(malnSet, comp->seq, comp->chromStart, comp->chromEnd, mafTreeLocAll);
    for (int i = 0; i < stList_length(overComps); i++) {
        struct malnComp *overComp = stList_get(overComps, i);
        if ((!malnDeleteBlks_contains(delBlks, overComp->blk))
             && checkMultiParent(comp, overComp, discardTwoParents)) {
            malnDeleteBlks_flag(delBlks, overComp->blk);
        }
    }
    stList_destruct(overComps);
}

/* check a block for components that have multiple parents  */
static void checkBlkForMultiParents(struct malnSet *malnSet, struct malnBlk *blk, struct malnDeleteBlks *delBlks, bool discardTwoParents) {
    blk->done = true;
    for (struct malnComp *comp = blk->comps; (comp != NULL); comp = comp->next) {
        if (malnComp_getLoc(comp) != mafTreeLocRoot) {
            checkCompForMultiParents(malnSet, comp, delBlks, discardTwoParents);
        }
    }
}

/* check for regions with multiple parents, deleting the blocks if requested,
 * otherwise aborting */
void malnMultiParents_check(struct malnSet *malnSet, bool discardTwoParents) {
    struct malnDeleteBlks *delBlks = malnDeleteBlks_construct();
    stSortedSetIterator *iter = malnSet_getBlocks(malnSet);
    struct malnBlk *blk;
    while ((blk = stSortedSet_getNext(iter)) != NULL) {
        if (!malnDeleteBlks_contains(delBlks, blk)) {
            checkBlkForMultiParents(malnSet, blk, delBlks, discardTwoParents);
        }
    }
    stSortedSet_destructIterator(iter);
    malnDeleteBlks_destruct(delBlks);
    malnSet_clearDone(malnSet);
}

