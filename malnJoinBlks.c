#include "malnJoinBlks.h"
#include "malnBlk.h"
#include "malnBlkCursor.h"
#include "malnComp.h"
#include "malnSet.h"
#include "mafTree.h"
#include "genome.h"
#include "common.h"
#include <stdbool.h>
#include <unistd.h>

static bool malnJoinBlksDebug = false;  // FIXME: tmp, make a flag

/*
 * join trees
 */
static mafTree *joinTrees(struct malnBlkCursor *blkCursor1, struct malnBlkCursor *blkCursor2) {
    struct malnComp *comp1 = blkCursor1->rows[0].comp;
    struct malnComp *comp2 = blkCursor2->rows[0].comp;
    return mafTree_join(blkCursor1->blk->mTree, comp1->seq->orgSeqName, comp1->chromStart, comp1->chromEnd,
                        blkCursor2->blk->mTree, comp2->seq->orgSeqName, comp2->chromStart, comp2->chromEnd);
}

/* create new reference components from the two being joined */
static struct malnComp *createJoinedRefComp(struct malnBlk *blkJoined, struct malnComp *refComp1, struct malnComp *refComp2) {
    int start = min(refComp1->start, refComp2->start);
    struct malnComp *comp = malnComp_construct(refComp1->seq, refComp1->strand, start, start, "");
    malnBlk_addComp(blkJoined, comp);
    return comp;
}

/* add non-reference components to the joined blk and destination component
 * array */
static void addCompsToJoined(struct malnBlk *blkJoined, struct malnBlkCursor *blkCursor, struct malnComp **destComps) {
    // skip rows[0], which is the reference
    for (int i = 1; i < blkCursor->numRows; i++) {
        destComps[i] = malnComp_construct(blkCursor->rows[i].comp->seq, blkCursor->rows[i].comp->strand, blkCursor->rows[i].comp->start, blkCursor->rows[i].comp->start, "");
        malnBlk_addComp(blkJoined, destComps[i]);
    }
}

/* assert new reference component covers the entire range */
static void assertJoinedRefComp(struct malnComp *destRefComp, struct malnComp *refComp1, struct malnComp *refComp2) {
    assert(destRefComp->start == min(refComp1->start, refComp2->start));
    assert(destRefComp->end == max(refComp1->end, refComp2->end));
}

/* assert new non-reference components cover entire range  */
static void assertJoinedComps(struct malnBlkCursor *blkCursor, struct malnComp **destComps) {
    // skip rows[0], which is the reference
    for (int i = 1; i < blkCursor->numRows; i++) {
        assert(destComps[i]->start == blkCursor->rows[i].comp->start);
        assert(destComps[i]->end == blkCursor->rows[i].comp->end);
    }
}

/* copy columns outside of the common reference sequence region to the joined maf */
static void copyUnsharedRefColumns(struct malnBlk *blkJoined, struct malnComp **destComps, struct malnBlkCursor *blkCursor, int alnStart, int alnEnd) {
    malnBlk_assert(blkJoined);  // FXIME
    malnBlk_assert(blkCursor->blk);  // FXIME
    if (false) {
        malnBlk_dump(blkJoined, "blkJoined", stderr);  // FXIME
        malnBlk_dump(blkCursor->blk, "blkCursor", stderr);  // FXIME
    }
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
    malnBlk_assert(blkJoined); // FIXME debugging
}

/* join two blocks using their specified reference components.  Optionally return
 * resulting join component. */
struct malnBlk *malnJoinBlks(struct malnComp *refComp1, struct malnComp *refComp2, struct malnComp **joinedCompRet) {
    assert(malnComp_overlap(refComp1, refComp2));
    struct malnBlk *blk1 = refComp1->blk;
    struct malnBlk *blk2 = refComp2->blk;
    if (malnJoinBlksDebug) { // FIXME: tmp
        malnBlk_dump(blk1, "inBlk1", stderr);
        malnBlk_dump(blk2, "inBlk2", stderr);
    }

    // reverse complement if needed
    struct malnBlk *freeBlk = NULL;
    if (refComp1->strand != refComp2->strand) {
        blk2 = freeBlk = malnBlk_reverseComplement(blk2);
        refComp2 = malnBlk_findCompByChromRange(blk2, refComp2->seq, refComp2->chromStart, refComp2->chromEnd);
        assert(refComp2 != NULL);
        if (malnJoinBlksDebug) { // FIXME: tmp
            malnBlk_dump(blk2, "inBlk2rc", stderr);
        }
    }

    // set up for the hard work destComps arrays parallel the rows in the
    // blkCursors
    struct malnBlkCursor *blkCursor1 = malnBlkCursor_construct(blk1, refComp1);
    struct malnBlkCursor *blkCursor2 = malnBlkCursor_construct(blk2, refComp2);
    struct malnComp *destComps1[blkCursor1->numRows];
    struct malnComp *destComps2[blkCursor2->numRows];
    struct malnBlk *blkJoined = malnBlk_construct(joinTrees(blkCursor1, blkCursor2));
    destComps1[0] = destComps2[0] = createJoinedRefComp(blkJoined, refComp1, refComp2);
    addCompsToJoined(blkJoined, blkCursor1, destComps1);
    addCompsToJoined(blkJoined, blkCursor2, destComps2);
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

    if (malnJoinBlksDebug) {
        fprintf(stderr, "refCommon: %s:%d-%d\n", refComp1->seq->orgSeqName, refCommonStart, refCommonEnd);
        fprintf(stderr, "align1Commom:: %d-%d\n", aln1CommonStart, aln1CommonEnd);
        fprintf(stderr, "align2Commom:: %d-%d\n", aln2CommonStart, aln2CommonEnd);
    }


    // before common start
    copyUnsharedRefColumns(blkJoined, destComps1, blkCursor1, 0, aln1CommonStart);
    if (malnJoinBlksDebug) {malnBlk_dump(blkJoined, "joined@1", stderr);}
    copyUnsharedRefColumns(blkJoined, destComps2, blkCursor2, 0, aln2CommonStart);
    if (malnJoinBlksDebug) {malnBlk_dump(blkJoined, "joined@2", stderr);}

    // common
    assert(destComps1[0]->end == refCommonStart);
    joinSharedRefColumns(blkJoined, destComps1, blkCursor1, aln1CommonEnd, destComps2, blkCursor2, aln2CommonEnd);
    assert(destComps1[0]->end == refCommonEnd);
    if (malnJoinBlksDebug) {malnBlk_dump(blkJoined, "joined@3", stderr);}

    // after common end
    copyUnsharedRefColumns(blkJoined, destComps1, blkCursor1, aln1CommonEnd, blk1->alnWidth);
    copyUnsharedRefColumns(blkJoined, destComps2, blkCursor2, aln2CommonEnd, blk2->alnWidth);
    
    assertJoinedRefComp(destComps1[0], refComp1, refComp2);
    assertJoinedComps(blkCursor1, destComps1);
    assertJoinedComps(blkCursor2, destComps2);
    malnBlkCursor_destruct(blkCursor1);
    malnBlkCursor_destruct(blkCursor2);
    blk1->done = blk2->done = true;  // must do before destruct, as blk1 might be freeBlk
    malnBlk_destruct(freeBlk);
    malnBlk_finish(blkJoined);
    malnBlk_assert(blkJoined);
    return blkJoined;
}
