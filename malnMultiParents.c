#include "malnMultiParents.h"
#include "malnSet.h"
#include "malnBlk.h"
#include "malnBlkMap.h"
#include "malnComp.h"
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
        fprintf(stderr, "%s\n", msg);
        malnBlk_dump(comp->blk, "compBlk",  stderr);
        malnBlk_dump(overComp->blk, "overCompBlk",  stderr);
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
static void checkCompForMultiParents(struct malnSet *malnSet, struct malnComp *comp, bool discardTwoParents, struct malnBlkMap *doneBlks) {
    stList *overComps = malnSet_getOverlappingPendingComps(malnSet, comp->seq, comp->chromStart, comp->chromEnd, mafTreeLocAll, doneBlks);
    for (int i = 0; i < stList_length(overComps); i++) {
        struct malnComp *overComp = stList_get(overComps, i);
        if (checkMultiParent(comp, overComp, discardTwoParents)) {
            malnBlk_markOrDelete(overComp->blk);
        }
    }
    stList_destruct(overComps);
}

/* check a block for components that have multiple parents  */
static void checkBlkForMultiParents(struct malnSet *malnSet, struct malnBlk *blk, bool discardTwoParents, struct malnBlkMap *doneBlks) {
    malnBlkMap_add(doneBlks, blk);
    for (struct malnComp *comp = blk->comps; (comp != NULL); comp = comp->next) {
        if (malnComp_getLoc(comp) != mafTreeLocRoot) {
            checkCompForMultiParents(malnSet, comp, discardTwoParents, doneBlks);
        }
    }
}

/* check for regions with multiple parents, deleting the blocks if requested,
 * otherwise aborting.  If deletingMap is not, skip entries in it.  */
void malnMultiParents_check(struct malnSet *malnSet, bool discardTwoParents) {
    struct malnBlkMap *doneBlks = malnBlkMap_construct();
    struct malnBlkMapIterator *iter = malnSet_getBlocks(malnSet);
    struct malnBlk *blk;
    while ((blk = malnBlkMapIterator_getNext(iter)) != NULL) {
        if (!blk->deleted) {
            checkBlkForMultiParents(malnSet, blk, discardTwoParents, doneBlks);
        }
    }
    malnSet_deleteDying(malnSet);
    malnBlkMapIterator_destruct(iter);
    malnBlkMap_destruct(doneBlks);
}

