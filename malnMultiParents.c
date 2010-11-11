#include "malnMultiParents.h"
#include "malnSet.h"
#include "malnBlk.h"
#include "malnBlkSet.h"
#include "malnComp.h"
#include "genome.h"
#include "sonLibList.h"
#include "stSafeC.h"

/* open file for reporting problems, writing header */
FILE *malnMultiParents_openResolveDropLog(char *multiParentDroppedFile) {
    static char *probHeader = "maf\tdroppedSeq\tdroppedStart\tdroppedEnd\n";
    FILE *dropLogFh = mustOpen(multiParentDroppedFile, "w");
    fputs(probHeader, dropLogFh);
    return dropLogFh;
}

/*
 * find overlapping coordinates of two components
 */
static int overlappingCoords(struct malnComp *comp1, struct malnComp *comp2, int *endRet) {
    assert(comp1->seq == comp2->seq);
    int start = max(comp1->chromStart, comp2->chromStart);
    int end = min(comp1->chromEnd, comp2->chromEnd);
    assert(start < end);
    *endRet = end;
    return start;
}

/* pick component to shorten.  returns 1 or 2. FIXME: make smarter,
 * probably want to remove the columns with the fewest aligned bases,
 * not just the shortest number of aligned columns */
static int pickCompToShorten(struct malnComp *comp1, int aln1Start, int aln1End, struct malnComp *comp2, int aln2Start, int aln2End) {
    if ((aln1End - aln1Start) < (aln2End - aln2Start)) {
        return 1;
    } else {
        return 2;
    }
}

/* shorten a component to eliminate common region */
static void shortenComp(struct malnComp *comp, int cutChromStart, int cutChromEnd) {
    // need to cut to one end, as middle cut isn't support (and maybe not desired)
    // leave the most sequence
    int chromStart, chromEnd;
    if ((cutChromStart == comp->chromStart) && (comp->chromEnd == cutChromEnd)) {
        // covers whole comp
        chromStart = comp->chromStart;
        chromEnd = comp->chromEnd;
    } else if ((cutChromStart - comp->chromStart) > (comp->chromEnd - cutChromEnd)) {
        // keep region at start
        chromStart = comp->chromStart;
        chromEnd = cutChromStart;
    } else {
        // keep region at end
        chromStart = cutChromEnd;
        chromEnd = comp->chromEnd;
    }
    malnBlk_shortenComp(comp->blk, comp, chromStart, chromEnd);
}

/* resolve conflicting parents for a component of a block.*/
static void resolveBlkOverComp(struct malnComp *comp1, struct malnComp *comp2, FILE *dropLogFh) {
    int commonChromStart, commonChromEnd;
    commonChromStart = overlappingCoords(comp1, comp2, &commonChromEnd);
    int aln1Start, aln1End, aln2Start, aln2End;
    malnComp_seqChromRangeToAlnRange(comp1, commonChromStart, commonChromEnd, &aln1Start, &aln1End);
    malnComp_seqChromRangeToAlnRange(comp2, commonChromStart, commonChromEnd, &aln2Start, &aln2End);

    if (dropLogFh != NULL) {
        fprintf(dropLogFh, "%s\t%s\t%d\t%d\n", malnSet_getMafFileName(comp1->blk->malnSet) ,comp1->seq->orgSeqName, commonChromStart, commonChromEnd);
    } else {
        fprintf(stderr, "Warning: multiple parents for %s:%d-%d, dropping one\n", comp1->seq->orgSeqName, commonChromStart, commonChromEnd);
    }

    if (pickCompToShorten(comp1, aln1Start, aln1End, comp2, aln2Start, aln2End) == 1) {
        shortenComp(comp1, commonChromStart, commonChromEnd);
    } else {
        shortenComp(comp2, commonChromStart, commonChromEnd);
    }
}

/* Resolve conflicting parents for a component of a block. Since comp
 * can be completely removed, can only update component and must start block
 * over.  Return true if something was done. */
static bool resolveBlkComp(struct malnSet *malnSet, struct malnBlk *blk, struct malnComp *comp, FILE *dropLogFh) {
    bool didSomething = false;
    stList *overComps = malnSet_getOverlappingPendingComps(malnSet, comp->seq, comp->chromStart, comp->chromEnd, mafTreeLocAll, NULL);
    for (int i = 0; (i < stList_length(overComps)) && (!didSomething) ; i++) {
        struct malnComp *overComp = stList_get(overComps, i);
        if ((overComp->blk != comp->blk) && (!overComp->blk->deleted) && (malnComp_getLoc(overComp) != mafTreeLocRoot)) {
            resolveBlkOverComp(comp, overComp, dropLogFh);
            didSomething = true;
        }
    }
    stList_destruct(overComps);
    return didSomething;
}

/* resolve conflicting parents on for block. */
static int resolveBlk(struct malnSet *malnSet, struct malnBlk *blk, FILE *dropLogFh) {
    // must start over each time something is resolved, since the block could
    // be updated.
    bool didSomething;
    int numResolved = 0;
    do {
        didSomething = false;
        for (struct malnComp *comp = blk->comps; comp != NULL; comp = comp->next) {
            if (malnComp_getLoc(comp) != mafTreeLocRoot) {
                didSomething = resolveBlkComp(malnSet, blk, comp, dropLogFh);
                if (didSomething) {
                    numResolved++;
                    break;  // start again
                }
            }
        }
    } while (didSomething);
    return numResolved;
}

/* Resolve conflicts where blocks imply sequences with multiple parents,
 * editing blocks.  If dropLogFh is not NULL, write a tab file describing the
 * problems */
int malnMultiParents_resolve(struct malnSet *malnSet, FILE *dropLogFh) {
    // this can be done in one pass, since the blocks and components are modified
    // in place instead of creating new ones
    int numResolved = 0;
    struct malnBlkSetIterator *iter = malnSet_getBlocks(malnSet);
    struct malnBlk *blk;
    while ((blk = malnBlkSetIterator_getNext(iter)) != NULL) {
        numResolved += resolveBlk(malnSet, blk, dropLogFh);
    }
    malnBlkSetIterator_destruct(iter);
    malnSet_deleteDying(malnSet);
    return numResolved;
}

/* report a multiple parent, either to stderr or by aborting */
static void reportMultiParent(struct malnComp *comp, struct malnComp *overComp) {
    char *msg = stSafeCDynFmt("multiple parents detected in components %s:%d-%d (%c) and %s:%d-%d (%c)",
                              comp->seq->orgSeqName, comp->start, comp->end, comp->strand,
                              overComp->seq->orgSeqName, overComp->start, overComp->end, overComp->strand);
    fprintf(stderr, "%s\n", msg);
    malnBlk_dump(comp->blk, stderr, "compBlk");
    malnBlk_dump(overComp->blk, stderr, "overCompBlk");
    errAbort("Error: %s", msg);
    stSafeCFree(msg);
}

/* compare a component to other components */
static void checkCompForMultiParents(struct malnSet *malnSet, struct malnComp *comp) {
    stList *overComps = malnSet_getOverlappingPendingComps(malnSet, comp->seq, comp->chromStart, comp->chromEnd, mafTreeLocAll, NULL);
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

