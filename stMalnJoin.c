#include "stMalnJoin.h"
#include "stMalnSet.h"
#include "stMalnBlk.h"
#include "stMalnBlkCursor.h"
#include "stMalnComp.h"
#include "stMafTree.h"
#include "genome.h"
#include "common.h"
#include <stdbool.h>
#include <unistd.h>

/* add components to the joined blk and destination component array, starting an row index,
 * to allow skipping reference */
static void addCompToJoined(struct stMalnBlk *blkJoined, int iStart, struct stMalnBlkCursor *blkCursor, struct stMalnComp **destComps) {
    for (int i = iStart; i < blkCursor->numRows; i++) {
        struct stMalnComp *comp = destComps[i] = stMalnComp_construct(blkCursor->rows[i].comp->seq, blkCursor->rows[i].comp->strand, blkCursor->rows[i].comp->start, blkCursor->rows[i].comp->start, "");
        stMalnBlk_addComp(blkJoined, comp);
    }
}

/*
 * join trees
 */
static stMafTree *joinTrees(struct stMalnBlkCursor *blkCursor1, struct stMalnBlkCursor *blkCursor2) {
    struct stMalnComp *comp1 = blkCursor1->rows[0].comp;
    struct stMalnComp *comp2 = blkCursor2->rows[0].comp;
    return stMafTree_join(blkCursor1->blk->mTree, comp1->seq->orgSeqName, comp1->chromStart, comp1->chromEnd,
                          blkCursor2->blk->mTree, comp2->seq->orgSeqName, comp2->chromStart, comp2->chromEnd);
}

/* initialize joined blk and an two array of destination components
 * corresponding to order in the two source block cursors, with the reference
 * sequence first in each */
static struct stMalnBlk *initJoinedBlk(struct stMalnBlkCursor *blkCursor1, struct stMalnBlkCursor *blkCursor2, struct stMalnComp ***destComps1Ret, struct stMalnComp ***destComps2Ret) {
    struct stMalnBlk *blkJoined = stMalnBlk_construct(joinTrees(blkCursor1, blkCursor2));
    struct stMalnComp **destComps1 = *destComps1Ret = needMem(blkCursor1->numRows*sizeof(struct stMalnComp*));
    struct stMalnComp **destComps2 = *destComps2Ret = needMem(blkCursor2->numRows*sizeof(struct stMalnComp*));

    // add components, reference is done only once
    addCompToJoined(blkJoined, 0, blkCursor1, destComps1);
    destComps2[0] = destComps1[0];
    addCompToJoined(blkJoined, 1, blkCursor2, destComps2);
    return blkJoined;
}

/* copy columns outside of the common reference sequence region to the joined maf */
static void copyUnsharedRefColumns(struct stMalnBlk *blkJoined, struct stMalnComp **destComps, struct stMalnBlkCursor *blkCursor, int alnStart, int alnEnd) {
    for (int i = 0; i < blkCursor->numRows; i++) {
        stMalnComp_appendCompAln(destComps[i], blkCursor->rows[i].comp, alnStart, alnEnd);
    }
    stMalnBlkCursor_setAlignCol(blkCursor, alnEnd);
    blkJoined->alnWidth += (alnEnd - alnStart);  // FIXME should have append methods
    stMalnBlk_pad(blkJoined);
    stMalnBlk_assert(blkJoined);
}

/* is reference aligned? */
static bool isRefAligned(struct stMalnBlkCursor *blkCursor) {
    return stMalnCompCursor_isAligned(&(blkCursor->rows[0]));
}

/* copy column from one source, optionally skipping reference */
static void copyColumn(struct stMalnComp **destComps, struct stMalnBlkCursor *blkCursor, bool skipRef) {
    for (int i = (skipRef ? 1 : 0); i < blkCursor->numRows; i++) {
        stMalnComp_appendCompAln(destComps[i], blkCursor->rows[i].comp, blkCursor->alnIdx, blkCursor->alnIdx+1);
    }
    stMalnBlkCursor_incr(blkCursor);
}

/* copy contiguous shared alignment columns to join blk */
static void copySharedRefColumns(struct stMalnBlk *blkJoined,
                                 struct stMalnComp **destComps1, struct stMalnBlkCursor *blkCursor1, int aln1CommonEnd,
                                 struct stMalnComp **destComps2, struct stMalnBlkCursor *blkCursor2, int aln2CommonEnd) {
    assert(isRefAligned(blkCursor1));
    assert(isRefAligned(blkCursor2));
    while (isRefAligned(blkCursor1) && isRefAligned(blkCursor2) && (blkCursor1->alnIdx < aln1CommonEnd) && (blkCursor2->alnIdx < aln2CommonEnd)) {
        copyColumn(destComps1, blkCursor1, false);
        copyColumn(destComps2, blkCursor2, true);
        blkJoined->alnWidth++;  // FIXME should have append methods
        stMalnBlk_assert(blkJoined);
    }
}

/* copy contiguous unaligned-to-reference columns to join blk */
static void copyUnalignedSharedColumns(struct stMalnBlk *blkJoined, struct stMalnComp **destComps, struct stMalnBlkCursor *blkCursor, int alnCommonEnd) {
    while ((!isRefAligned(blkCursor)) && (blkCursor->alnIdx < alnCommonEnd)) {
        copyColumn(destComps, blkCursor, false);
        blkJoined->alnWidth++;  // FIXME should have append methods
    }
    stMalnBlk_pad(blkJoined);
    stMalnBlk_assert(blkJoined);
}

/* join columns based on shared reference sequence regions */
static void joinSharedRefColumns(struct stMalnBlk *blkJoined,
                                 struct stMalnComp **destComps1, struct stMalnBlkCursor *blkCursor1, int aln1CommonEnd,
                                 struct stMalnComp **destComps2, struct stMalnBlkCursor *blkCursor2, int aln2CommonEnd) {
    assert(blkCursor1->rows[0].pos == blkCursor2->rows[0].pos);
    while ((blkCursor1->alnIdx < aln1CommonEnd) && (blkCursor1->alnIdx < aln1CommonEnd)) {
        copySharedRefColumns(blkJoined, destComps1, blkCursor1, aln1CommonEnd, destComps2, blkCursor2, aln2CommonEnd);
        copyUnalignedSharedColumns(blkJoined, destComps1, blkCursor1, aln1CommonEnd);
        copyUnalignedSharedColumns(blkJoined, destComps2, blkCursor2, aln2CommonEnd);
    }
    assert(blkCursor1->alnIdx == aln1CommonEnd);
    assert(blkCursor2->alnIdx == aln2CommonEnd);
    stMalnBlk_assert(blkJoined);
}

/* join two blocks using the specified reference components. */
static struct stMalnBlk *joinBlksByRef(struct stMalnBlk *blk1, struct stMalnComp *refComp1, struct stMalnBlk *blk2, struct stMalnComp *refComp2) {
    // reverse complement if needed
    struct stMalnBlk *freeBlk = NULL;
    if (refComp1->strand != refComp2->strand) {
        freeBlk = stMalnBlk_reverseComplement(blk2);
        blk2 = freeBlk;
        refComp2 = stMalnBlk_findCompByChromRange(blk2, refComp2->seq, refComp2->chromStart, refComp2->chromEnd);
        assert(refComp2 != NULL);
    }

    // set up for the hard work
    struct stMalnBlkCursor *blkCursor1 = stMalnBlkCursor_construct(blk1, refComp1);
    struct stMalnBlkCursor *blkCursor2 = stMalnBlkCursor_construct(blk2, refComp2);
    struct stMalnComp **destComps1, **destComps2;
    struct stMalnBlk *blkJoined = initJoinedBlk(blkCursor1, blkCursor2, &destComps1, &destComps2);
    int refCommonStart = max(refComp1->start, refComp2->start);
    int refCommonEnd = min(refComp1->end, refComp2->end);
    int aln1CommonStart, aln1CommonEnd, aln2CommonStart, aln2CommonEnd;
    if (!(stMalnComp_seqRangeToAlnRange(refComp1, refCommonStart, refCommonEnd, &aln1CommonStart, &aln1CommonEnd)
          && stMalnComp_seqRangeToAlnRange(refComp2, refCommonStart, refCommonEnd, &aln2CommonStart, &aln2CommonEnd))) {
        errAbort("BUG: failure to get alignment ranges for common reference sequence range");
    }
    
    // before common start
    copyUnsharedRefColumns(blkJoined, destComps1, blkCursor1, 0, aln1CommonStart);
    copyUnsharedRefColumns(blkJoined, destComps2, blkCursor2, 0, aln2CommonStart);

    // common
    joinSharedRefColumns(blkJoined, destComps1, blkCursor1, aln1CommonEnd, destComps2, blkCursor2, aln2CommonEnd);

    // after common end
    copyUnsharedRefColumns(blkJoined, destComps1, blkCursor1, aln1CommonEnd, blk1->alnWidth);
    copyUnsharedRefColumns(blkJoined, destComps2, blkCursor2, aln2CommonEnd, blk2->alnWidth);

    freeMem(destComps1);
    freeMem(destComps2);
    stMalnBlkCursor_destruct(blkCursor1);
    stMalnBlkCursor_destruct(blkCursor2);
    stMalnBlk_destruct(freeBlk);
    blk1->done = blk2->done = true;
    stMalnBlk_sortComps(blkJoined);
    stMalnBlk_assert(blkJoined);
    return blkJoined;
}

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
        struct stMalnBlk *blkJoined = joinBlksByRef(blk1, refComp1, blk2, refComp2);
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
