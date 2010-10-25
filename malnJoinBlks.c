#include "malnJoinBlks.h"
#include "malnBlk.h"
#include "malnBlkCursor.h"
#include "malnComp.h"
#include "malnSet.h"
#include "mafTree.h"
#include "malnCompCompMap.h"
#include "genome.h"
#include "common.h"
#include <stdbool.h>
#include <unistd.h>

/* Object use to keep state */
struct malnJoinBlks {
    struct malnComp *guide1;         // guide components
    struct malnComp *guide2;
    struct malnBlk *blk1;            // input blocks
    struct malnBlk *blk2;
    struct malnBlk *inBlk2;          // original blk2, before reverse complement
    struct malnBlkCursor *cursor1;   // cursors into input block
    struct malnBlkCursor *cursor2;
    struct malnBlk *joined;          // joined block
    struct malnComp **dests1;        // arrays of source to destination components, first common
    struct malnComp **dests2;
    struct malnCompCompMap *srcDestCompMap; // map of source to destination components for tree join
    struct malnBlk *freeBlk;   // block to free if not NULL (due to reverse-complement);

    // coordinates of the two alignments
    int guideCommonStart;  // common guide sequence coordinates, zero length for adjacent
    int guideCommonEnd;
    int aln1CommonStart; // mapping of the guide coordinates to input alignments
    int aln1CommonEnd;
    int aln2CommonStart;
    int aln2CommonEnd;
};

/* compute alignment coordinates */
static void calcAlignmentCoords(struct malnJoinBlks *jb) {
    jb->guideCommonStart = max(jb->guide1->start, jb->guide2->start);
    jb->guideCommonEnd = min(jb->guide1->end, jb->guide2->end);
    if (jb->guideCommonStart < jb->guideCommonEnd) {
        // overlapping
        if (!(malnComp_seqRangeToAlnRange(jb->guide1, jb->guideCommonStart, jb->guideCommonEnd, &jb->aln1CommonStart, &jb->aln1CommonEnd)
              && malnComp_seqRangeToAlnRange(jb->guide2, jb->guideCommonStart, jb->guideCommonEnd, &jb->aln2CommonStart, &jb->aln2CommonEnd))) {
            errAbort("BUG: failure to get alignment ranges for common guide sequence range");
        }
    } else {
        // adjacent
        jb->aln1CommonStart = jb->aln1CommonEnd = jb->cursor1->alnWidth;
        jb->aln2CommonStart = jb->aln2CommonEnd = jb->cursor2->alnWidth;
    }
}

/* create new guide components from the two being joined */
static struct malnComp *createJoinedGuideComp(struct malnJoinBlks *jb) {
    int start = min(jb->guide1->start, jb->guide2->start);
    struct malnComp *comp = malnComp_construct(jb->guide1->seq, jb->guide1->strand, start, start, "");
    malnBlk_addComp(jb->joined, comp);
    malnCompCompMap_add(jb->srcDestCompMap, jb->guide1, comp);
    malnCompCompMap_add(jb->srcDestCompMap, jb->guide2, comp);
    return comp;
}

/* add non-guide components to the joined blk and destination component
 * array */
static void addCompsToJoined(struct malnJoinBlks *jb, struct malnBlkCursor *blkCursor, struct malnComp **destComps) {
    // skip rows[0], which is the guide
    for (int i = 1; i < blkCursor->numRows; i++) {
        destComps[i] = malnComp_construct(blkCursor->rows[i].comp->seq, blkCursor->rows[i].comp->strand, blkCursor->rows[i].comp->start, blkCursor->rows[i].comp->start, "");
        malnBlk_addComp(jb->joined, destComps[i]);
        malnCompCompMap_add(jb->srcDestCompMap, blkCursor->rows[i].comp, destComps[i]);
    }
}

/* construct malnJoinBlks state object for the join */
static struct malnJoinBlks *malnJoinBlks_construct(struct malnComp *guideComp1, struct malnComp *guideComp2) {
    struct malnJoinBlks *jb;
    AllocVar(jb);
    jb->guide1 = guideComp1;
    jb->guide2 = guideComp2;
    jb->blk1 = guideComp1->blk;
    jb->blk2 = jb->inBlk2 = guideComp2->blk;

    // reverse complement if needed (must do before making cursors)
    if (jb->guide1->strand != jb->guide2->strand) {
        jb->blk2 = jb->freeBlk = malnBlk_reverseComplement(jb->blk2);
        jb->guide2 = malnBlk_findCompByChromRange(jb->blk2, jb->guide2->seq, jb->guide2->chromStart, jb->guide2->chromEnd);
        assert(jb->guide2 != NULL);
    }

    jb->cursor1 = malnBlkCursor_construct(jb->blk1, jb->guide1, NULL);
    jb->cursor2 = malnBlkCursor_construct(jb->blk2, jb->guide2, NULL);
    jb->joined = malnBlk_construct();
    jb->dests1 = needMem(jb->cursor1->numRows * sizeof(struct malnComp *));
    jb->dests2 = needMem(jb->cursor2->numRows * sizeof(struct malnComp *));
    jb->srcDestCompMap = malnCompCompMap_construct();
    jb->freeBlk = NULL;

    jb->dests1[0] = jb->dests2[0] = createJoinedGuideComp(jb);
    addCompsToJoined(jb, jb->cursor1, jb->dests1);
    addCompsToJoined(jb, jb->cursor2, jb->dests2);
    return jb;
}

/* destruct state object */
static void malnJoinBlks_destruct(struct malnJoinBlks *jb) {
    malnCompCompMap_destruct(jb->srcDestCompMap);
    malnBlkCursor_destruct(jb->cursor1);
    malnBlkCursor_destruct(jb->cursor2);
    malnBlk_destruct(jb->freeBlk);
    freeMem(jb->dests1);
    freeMem(jb->dests2);
    freeMem(jb);
}

/*
 * join trees
 */
static mafTree *joinTrees(struct malnBlkCursor *blkCursor1, struct malnBlkCursor *blkCursor2, struct malnCompCompMap *srcDestCompMap) {
    struct malnComp *comp1 = blkCursor1->rows[0].comp;
    struct malnComp *comp2 = blkCursor2->rows[0].comp;
    return mafTree_join(blkCursor1->blk->mTree, comp1, blkCursor2->blk->mTree, comp2, srcDestCompMap);
}

/* assert new guide component covers the entire range */
static void assertJoinedGuideComp(struct malnComp *destGuideComp, struct malnComp *guideComp1, struct malnComp *guideComp2) {
    assert(destGuideComp->start == min(guideComp1->start, guideComp2->start));
    assert(destGuideComp->end == max(guideComp1->end, guideComp2->end));
}

/* assert new non-guide components cover entire range  */
static void assertJoinedComps(struct malnBlkCursor *blkCursor, struct malnComp **destComps) {
    // skip rows[0], which is the guide
    for (int i = 1; i < blkCursor->numRows; i++) {
        assert(destComps[i]->start == blkCursor->rows[i].comp->start);
        assert(destComps[i]->end == blkCursor->rows[i].comp->end);
    }
}

/* copy columns outside of the common guide sequence region to the joined maf */
static void copyUnsharedGuideColumns(struct malnJoinBlks *jb, struct malnComp **destComps, struct malnBlkCursor *blkCursor, int alnStart, int alnEnd) {
  for (int i = 0; i < blkCursor->numRows; i++) {
        malnComp_appendCompAln(destComps[i], blkCursor->rows[i].comp, alnStart, alnEnd);
    }
    malnBlkCursor_setAlignCol(blkCursor, alnEnd);
    jb->joined->alnWidth += (alnEnd - alnStart);  // FIXME should have append methods
    malnBlk_pad(jb->joined);
}

/* is guide aligned? */
static bool isGuideAligned(struct malnBlkCursor *blkCursor) {
    return malnCompCursor_isAligned(&(blkCursor->rows[0]));
}

/* copy column from one source, optionally skipping guide */
static void copyColumn(struct malnComp **destComps, struct malnBlkCursor *blkCursor, bool skipGuide) {
    for (int i = (skipGuide ? 1 : 0); i < blkCursor->numRows; i++) {
        malnComp_appendCompAln(destComps[i], blkCursor->rows[i].comp, blkCursor->alnIdx, blkCursor->alnIdx+1);
    }
    malnBlkCursor_incr(blkCursor);
}

/* copy contiguous shared alignment columns to join blk */
static void copySharedGuideColumns(struct malnJoinBlks *jb) {
    assert(isGuideAligned(jb->cursor1));
    assert(isGuideAligned(jb->cursor2));
    while (isGuideAligned(jb->cursor1) && isGuideAligned(jb->cursor2) && (jb->cursor1->alnIdx < jb->aln1CommonEnd) && (jb->cursor2->alnIdx < jb->aln2CommonEnd)) {
        copyColumn(jb->dests1, jb->cursor1, false);
        copyColumn(jb->dests2, jb->cursor2, true);
        jb->joined->alnWidth++;
    }
}

/* copy contiguous unaligned-to-guide columns to join blk */
static void copyUnalignedSharedColumns(struct malnBlk *blkJoined, struct malnComp **destComps, struct malnBlkCursor *blkCursor, int alnCommonEnd) {
    while ((!isGuideAligned(blkCursor)) && (blkCursor->alnIdx < alnCommonEnd)) {
        copyColumn(destComps, blkCursor, false);
        blkJoined->alnWidth++;
    }
    malnBlk_pad(blkJoined);
}

/* join columns based on shared guide sequence regions */
static void joinSharedGuideColumns(struct malnJoinBlks *jb) {
    assert(jb->cursor1->rows[0].pos == jb->cursor2->rows[0].pos);
    while ((jb->cursor1->alnIdx < jb->aln1CommonEnd) && (jb->cursor1->alnIdx < jb->aln1CommonEnd)) {
        copySharedGuideColumns(jb);
        copyUnalignedSharedColumns(jb->joined, jb->dests1, jb->cursor1, jb->aln1CommonEnd);
        copyUnalignedSharedColumns(jb->joined, jb->dests2, jb->cursor2, jb->aln2CommonEnd);
    }
    assert(jb->cursor1->alnIdx == jb->aln1CommonEnd);
    assert(jb->cursor2->alnIdx == jb->aln2CommonEnd);
}

/* join two blocks using their specified guide components.  Optionally return
 * resulting join component. */
struct malnBlk *malnJoinBlks(struct malnComp *guideComp1, struct malnComp *guideComp2, struct malnComp **joinedCompRet) {
    // FIXME: is joinedCompRet needed??
    assert(malnComp_overlapAdjacent(guideComp1, guideComp2));
    struct malnJoinBlks *jb = malnJoinBlks_construct(guideComp1, guideComp2);
    if (joinedCompRet != NULL) {
        *joinedCompRet = jb->dests1[0];
    }
    calcAlignmentCoords(jb);

    // before common start
    copyUnsharedGuideColumns(jb, jb->dests1, jb->cursor1, 0, jb->aln1CommonStart);
    copyUnsharedGuideColumns(jb, jb->dests2, jb->cursor2, 0, jb->aln2CommonStart);

    // common
    if (jb->guideCommonStart < jb->guideCommonEnd) {
        assert(jb->dests1[0]->end == jb->guideCommonStart);
        joinSharedGuideColumns(jb);
        assert(jb->dests1[0]->end == jb->guideCommonEnd);
    }

    // after common end
    copyUnsharedGuideColumns(jb, jb->dests1, jb->cursor1, jb->aln1CommonEnd, jb->blk1->alnWidth);
    copyUnsharedGuideColumns(jb, jb->dests2, jb->cursor2, jb->aln2CommonEnd, jb->blk2->alnWidth);

    malnBlk_setTree(jb->joined, joinTrees(jb->cursor1, jb->cursor2, jb->srcDestCompMap));

    assertJoinedGuideComp(jb->dests1[0], jb->guide1, jb->guide2);
    assertJoinedComps(jb->cursor1, jb->dests1);
    assertJoinedComps(jb->cursor2, jb->dests2);
    malnBlk_finish(jb->joined);

    struct malnBlk *joined = jb->joined;
    malnJoinBlks_destruct(jb);
    return joined;
}
