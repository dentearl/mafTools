#include "malnJoinBlks.h"
#include "malnBlk.h"
#include "malnBlkCursor.h"
#include "malnComp.h"
#include "mafTree.h"
#include "genome.h"
#include "common.h"
#include <stdbool.h>
#include <unistd.h>

/*
 * join trees
 */
static mafTree *joinTrees(struct malnBlkCursor *blkCursor1, struct malnBlkCursor *blkCursor2) {
    struct malnComp *comp1 = blkCursor1->rows[0].comp;
    struct malnComp *comp2 = blkCursor2->rows[0].comp;
    return mafTree_join(blkCursor1->blk->mTree, comp1->seq->orgSeqName, comp1->chromStart, comp1->chromEnd,
                        blkCursor2->blk->mTree, comp2->seq->orgSeqName, comp2->chromStart, comp2->chromEnd);
}

/* add components to the joined blk and destination component array, starting an row index,
 * to allow skipping reference */
static void addCompToJoined(struct malnBlk *blkJoined, int iStart, struct malnBlkCursor *blkCursor, struct malnComp **destComps) {
    for (int i = iStart; i < blkCursor->numRows; i++) {
        struct malnComp *comp = destComps[i] = malnComp_construct(blkCursor->rows[i].comp->seq, blkCursor->rows[i].comp->strand, blkCursor->rows[i].comp->start, blkCursor->rows[i].comp->start, "");
        malnBlk_addComp(blkJoined, comp);
    }
}

/* initialize joined blk and an two array of destination components
 * corresponding to order in the two source block cursors, with the reference
 * sequence first in each */
static struct malnBlk *initJoinedBlk(struct malnBlkCursor *blkCursor1, struct malnBlkCursor *blkCursor2, struct malnComp ***destComps1Ret, struct malnComp ***destComps2Ret) {
    struct malnBlk *blkJoined = malnBlk_construct(joinTrees(blkCursor1, blkCursor2));
    struct malnComp **destComps1 = *destComps1Ret = needMem(blkCursor1->numRows*sizeof(struct malnComp*));
    struct malnComp **destComps2 = *destComps2Ret = needMem(blkCursor2->numRows*sizeof(struct malnComp*));

    // add components, reference is done only once
    addCompToJoined(blkJoined, 0, blkCursor1, destComps1);
    destComps2[0] = destComps1[0];
    addCompToJoined(blkJoined, 1, blkCursor2, destComps2);
    return blkJoined;
}

/* copy columns outside of the common reference sequence region to the joined maf */
static void copyUnsharedRefColumns(struct malnBlk *blkJoined, struct malnComp **destComps, struct malnBlkCursor *blkCursor, int alnStart, int alnEnd) {
    for (int i = 0; i < blkCursor->numRows; i++) {
        malnComp_appendCompAln(destComps[i], blkCursor->rows[i].comp, alnStart, alnEnd);
    }
    malnBlkCursor_setAlignCol(blkCursor, alnEnd);
    blkJoined->alnWidth += (alnEnd - alnStart);  // FIXME should have append methods
    malnBlk_pad(blkJoined);
    malnBlk_assert(blkJoined);
}

/* is reference aligned? */
static bool isRefAligned(struct malnBlkCursor *blkCursor) {
    return malnCompCursor_isAligned(&(blkCursor->rows[0]));
}

/* copy column from one source, optionally skipping reference */
static void copyColumn(struct malnComp **destComps, struct malnBlkCursor *blkCursor, bool skipRef) {
    for (int i = (skipRef ? 1 : 0); i < blkCursor->numRows; i++) {
        malnComp_appendCompAln(destComps[i], blkCursor->rows[i].comp, blkCursor->alnIdx, blkCursor->alnIdx+1);
    }
    malnBlkCursor_incr(blkCursor);
}

/* copy contiguous shared alignment columns to join blk */
static void copySharedRefColumns(struct malnBlk *blkJoined,
                                 struct malnComp **destComps1, struct malnBlkCursor *blkCursor1, int aln1CommonEnd,
                                 struct malnComp **destComps2, struct malnBlkCursor *blkCursor2, int aln2CommonEnd) {
    assert(isRefAligned(blkCursor1));
    assert(isRefAligned(blkCursor2));
    while (isRefAligned(blkCursor1) && isRefAligned(blkCursor2) && (blkCursor1->alnIdx < aln1CommonEnd) && (blkCursor2->alnIdx < aln2CommonEnd)) {
        copyColumn(destComps1, blkCursor1, false);
        copyColumn(destComps2, blkCursor2, true);
        blkJoined->alnWidth++;  // FIXME should have append methods
        malnBlk_assert(blkJoined);
    }
}

/* copy contiguous unaligned-to-reference columns to join blk */
static void copyUnalignedSharedColumns(struct malnBlk *blkJoined, struct malnComp **destComps, struct malnBlkCursor *blkCursor, int alnCommonEnd) {
    while ((!isRefAligned(blkCursor)) && (blkCursor->alnIdx < alnCommonEnd)) {
        copyColumn(destComps, blkCursor, false);
        blkJoined->alnWidth++;  // FIXME should have append methods
    }
    malnBlk_pad(blkJoined);
    malnBlk_assert(blkJoined);
}

/* join columns based on shared reference sequence regions */
static void joinSharedRefColumns(struct malnBlk *blkJoined,
                                 struct malnComp **destComps1, struct malnBlkCursor *blkCursor1, int aln1CommonEnd,
                                 struct malnComp **destComps2, struct malnBlkCursor *blkCursor2, int aln2CommonEnd) {
    assert(blkCursor1->rows[0].pos == blkCursor2->rows[0].pos);
    while ((blkCursor1->alnIdx < aln1CommonEnd) && (blkCursor1->alnIdx < aln1CommonEnd)) {
        copySharedRefColumns(blkJoined, destComps1, blkCursor1, aln1CommonEnd, destComps2, blkCursor2, aln2CommonEnd);
        copyUnalignedSharedColumns(blkJoined, destComps1, blkCursor1, aln1CommonEnd);
        copyUnalignedSharedColumns(blkJoined, destComps2, blkCursor2, aln2CommonEnd);
    }
    assert(blkCursor1->alnIdx == aln1CommonEnd);
    assert(blkCursor2->alnIdx == aln2CommonEnd);
    malnBlk_assert(blkJoined);
}

/* join two blocks using their specified reference components.  Optionally return
 * resulting join component. */
struct malnBlk *malnJoinBlks(struct malnComp *refComp1, struct malnComp *refComp2, struct malnComp **joinedCompRet) {
    assert(malnComp_overlap(refComp1, refComp2));
    struct malnBlk *blk1 = refComp1->blk;
    struct malnBlk *blk2 = refComp2->blk;

    // reverse complement if needed
    struct malnBlk *freeBlk = NULL;
    if (refComp1->strand != refComp2->strand) {
        freeBlk = malnBlk_reverseComplement(blk2);
        blk2 = freeBlk;
        refComp2 = malnBlk_findCompByChromRange(blk2, refComp2->seq, refComp2->chromStart, refComp2->chromEnd);
        assert(refComp2 != NULL);
    }

    // set up for the hard work
    struct malnBlkCursor *blkCursor1 = malnBlkCursor_construct(blk1, refComp1);
    struct malnBlkCursor *blkCursor2 = malnBlkCursor_construct(blk2, refComp2);
    struct malnComp **destComps1, **destComps2;
    struct malnBlk *blkJoined = initJoinedBlk(blkCursor1, blkCursor2, &destComps1, &destComps2);
    int refCommonStart = max(refComp1->start, refComp2->start);
    int refCommonEnd = min(refComp1->end, refComp2->end);
    int aln1CommonStart, aln1CommonEnd, aln2CommonStart, aln2CommonEnd;
    if (!(malnComp_seqRangeToAlnRange(refComp1, refCommonStart, refCommonEnd, &aln1CommonStart, &aln1CommonEnd)
          && malnComp_seqRangeToAlnRange(refComp2, refCommonStart, refCommonEnd, &aln2CommonStart, &aln2CommonEnd))) {
        errAbort("BUG: failure to get alignment ranges for common reference sequence range");
    }
    if (joinedCompRet != NULL) {
        *joinedCompRet = destComps1[0];
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
    malnBlkCursor_destruct(blkCursor1);
    malnBlkCursor_destruct(blkCursor2);
    malnBlk_destruct(freeBlk);
    blk1->done = blk2->done = true;
    malnBlk_sortComps(blkJoined);
    malnBlk_setLocAttr(blkJoined, refComp1->seq->genome);
    malnBlk_assert(blkJoined);
    return blkJoined;
}


