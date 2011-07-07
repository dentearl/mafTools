#include "malnMultiParents.h"
#include "malnSet.h"
#include "malnBlk.h"
#include "malnBlkSet.h"
#include "malnComp.h"
#include "genome.h"
#include "sonLibList.h"
#include "stSafeC.h"

/* report a multiple parent, either to stderr or by aborting */
static void reportMultiParent(struct malnComp *comp, struct malnComp *overComp) {
    char *msg = stSafeCDynFmt("multiple parents detected in components %s:%d-%d (%c) and %s:%d-%d (%c)",
                              comp->seq->orgSeqName, comp->start, comp->end, comp->strand,
                              overComp->seq->orgSeqName, overComp->start, overComp->end, overComp->strand);
    fprintf(stderr, "%s\n", msg);
    malnBlk_dump(comp->blk, stderr, malnBlk_dumpDefault, "compBlk");
    malnBlk_dump(overComp->blk, stderr, malnBlk_dumpDefault, "overCompBlk");
    errAbort("Error: %s", msg);
    stSafeCFree(msg);
}

/* compare a component to other components */
static void checkCompForMultiParents(struct malnSet *malnSet, struct malnComp *comp) {
    stList *overComps = malnSet_getOverlappingComps(malnSet, comp->seq, comp->chromStart, comp->chromEnd, mafTreeLocAll);
    for (int i = 0; i < stList_length(overComps); i++) {
        struct malnComp *overComp = stList_get(overComps, i);
        if ((comp->blk != overComp->blk) && (malnComp_getLoc(overComp) != mafTreeLocRoot)) {
            reportMultiParent(comp, overComp);
        }
    }
    stList_destruct(overComps);
}

/* check a block for components that have multiple parents  */
static void checkBlkForMultiParents(struct malnSet *malnSet, struct malnBlk *blk) {
    for (struct malnComp *comp = blk->comps; (comp != NULL); comp = comp->next) {
        if (malnComp_getLoc(comp) != mafTreeLocRoot) {
            checkCompForMultiParents(malnSet, comp);
        }
    }
}

/* check for regions with multiple parents, deleting the blocks if requested,
 * otherwise aborting. */
void malnMultiParents_check(struct malnSet *malnSet) {
    struct malnBlkSetIterator *iter = malnSet_getBlocks(malnSet);
    struct malnBlk *blk;
    while ((blk = malnBlkSetIterator_getNext(iter)) != NULL) {
        if (!blk->deleted) {
            checkBlkForMultiParents(malnSet, blk);
        }
    }
    malnBlkSetIterator_destruct(iter);
}

