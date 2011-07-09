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

static bool debug = false; // FIXME: tmp
//static unsigned debugDumpOpts = malnBlk_dumpNoSeqs ; // FIXME: tmp
static unsigned debugDumpOpts = malnBlk_dumpDefault ; // FIXME: tmp

/* Object use to keep state */
struct malnJoinBlks {
    struct malnComp *guide1;         // guide components, #1 start before or at #2
    struct malnComp *guide2;
    struct malnBlk *blk1;            // input blocks
    struct malnBlk *blk2;
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
    struct malnComp *comp = malnComp_construct(jb->guide1->seq, jb->guide1->strand, start, start, 0);
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
        destComps[i] = malnComp_construct(blkCursor->rows[i].comp->seq, blkCursor->rows[i].comp->strand, blkCursor->rows[i].comp->start, blkCursor->rows[i].comp->start, 0);
        malnBlk_addComp(jb->joined, destComps[i]);
        malnCompCompMap_add(jb->srcDestCompMap, blkCursor->rows[i].comp, destComps[i]);
    }
}

/* reverse complement one of the blocks give the guide component */
static struct malnComp *revComplementGuide(struct malnComp *guideComp, struct malnBlk **freeBlk) {
    *freeBlk = malnBlk_reverseComplement(guideComp->blk);
    guideComp = malnBlk_findCompByChromRange(*freeBlk, guideComp->seq, guideComp->chromStart, guideComp->chromEnd);
    assert(guideComp != NULL);
    return guideComp;
}

/* construct malnJoinBlks state object for the join */
static struct malnJoinBlks *malnJoinBlks_construct(struct malnComp *guideComp1, struct malnComp *guideComp2) {
    struct malnJoinBlks *jb;
    AllocVar(jb);
    // reverse complement a guide if needed (must do first), pick so they are
    // positive strand
    if (guideComp1->strand != guideComp2->strand) {
        if (guideComp1->strand == '-') {
            guideComp1 = revComplementGuide(guideComp1, &jb->freeBlk);
        } else {
            guideComp2 = revComplementGuide(guideComp2, &jb->freeBlk);
        }
    }

    // make sure guideComp1 starts before or at guideComp2 (although don't
    // remember why)
    if (guideComp1->start > guideComp2->start) {
        struct malnComp *hold = guideComp2;
        guideComp2 = guideComp1;
        guideComp1 = hold;
    }
    jb->guide1 = guideComp1;
    jb->guide2 = guideComp2;
    jb->blk1 = guideComp1->blk;
    jb->blk2 = guideComp2->blk;

    jb->cursor1 = malnBlkCursor_construct(jb->blk1, jb->guide1, NULL);
    malnBlkCursor_incr(jb->cursor1); // advance to start
    jb->cursor2 = malnBlkCursor_construct(jb->blk2, jb->guide2, NULL);
    malnBlkCursor_incr(jb->cursor2);
    jb->joined = malnBlk_construct();
    jb->dests1 = needMem(jb->cursor1->numRows * sizeof(struct malnComp *));
    jb->dests2 = needMem(jb->cursor2->numRows * sizeof(struct malnComp *));
    jb->srcDestCompMap = malnCompCompMap_construct();

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
    assert(alnStart == blkCursor->alnIdx);  // FIXME: alnStart not really needed
    for (int i = 0; i < blkCursor->numRows; i++) {
        malnComp_appendFromCursor(destComps[i], &(blkCursor->rows[i]), alnEnd);
    }
    malnBlkCursor_setAlignCol(blkCursor, alnEnd);
    malnBlk_pad(jb->joined);
}

#ifndef NDEBUG
/* is guide aligned? */
static bool isGuideAligned(struct malnBlkCursor *blkCursor) {
    return malnCompCursor_isAligned(&(blkCursor->rows[0]));
}
#endif

/* get end index of contiguous aligned or unaligned guide sequence columns, up to the specified end point */
static int getContiguousGuideStateEnd(struct malnBlkCursor *blkCursor, int alnEnd, bool wantAligned) {
    // FIXME: way too inefficient, should use range.
    struct malnCompCursor guideCursor = blkCursor->rows[0];
    while ((guideCursor.alnIdx < alnEnd) && (malnCompCursor_isAligned(&guideCursor) == wantAligned)) {
        if (!malnCompCursor_incr(&guideCursor)) {
            break;
        }
    }
    return guideCursor.alnIdx;
}

/* copy columns from one alignment to another, updating the cursor. */
static void copyColumns(struct malnComp **destComps, struct malnBlkCursor *blkCursor, int alnEnd, bool skipGuide) {
    if (skipGuide) {
        malnCompCursor_setAlignCol(&(blkCursor->rows[0]), alnEnd);
    }
    for (int i = (skipGuide ? 1 : 0); i < blkCursor->numRows; i++) {
        malnComp_appendFromCursor(destComps[i], &(blkCursor->rows[i]), alnEnd);
    }
    blkCursor->alnIdx = alnEnd;    // FIXME; make function
}

/* copy contiguous shared alignment columns to join blk */
static void copySharedGuideColumns(struct malnJoinBlks *jb) {
    assert(isGuideAligned(jb->cursor1));
    assert(isGuideAligned(jb->cursor2));
    // get minimum of the two that are aligned to the comm
    int aln1End = getContiguousGuideStateEnd(jb->cursor1, jb->aln1CommonEnd, true);
    int aln2End = getContiguousGuideStateEnd(jb->cursor2, jb->aln2CommonEnd, true);
    int numCols = min((aln1End - jb->cursor1->alnIdx), (aln2End - jb->cursor2->alnIdx));

    copyColumns(jb->dests1, jb->cursor1, jb->cursor1->alnIdx + numCols, false);
    copyColumns(jb->dests2, jb->cursor2, jb->cursor2->alnIdx + numCols, true);
}

/* copy contiguous unaligned-to-guide columns to join blk */
static void copyUnalignedSharedColumns(struct malnBlk *blkJoined, struct malnComp **destComps, struct malnBlkCursor *blkCursor, int alnEnd) {
    alnEnd = getContiguousGuideStateEnd(blkCursor, alnEnd, false);
    copyColumns(destComps, blkCursor, alnEnd, false);
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

/* join two blocks using their specified guide components. */
struct malnBlk *malnJoinBlks(struct malnComp *guideComp1, struct malnComp *guideComp2) {
    assert(malnComp_overlapAdjacent(guideComp1, guideComp2));
    struct malnJoinBlks *jb = malnJoinBlks_construct(guideComp1, guideComp2);
    calcAlignmentCoords(jb);
    if (debug) {
        malnBlk_dump(jb->blk1, stderr, debugDumpOpts, "inBlk1");
        malnBlk_dump(jb->blk2, stderr, debugDumpOpts, "inBlk2");
    }

    // before common start
    copyUnsharedGuideColumns(jb, jb->dests1, jb->cursor1, 0, jb->aln1CommonStart);
    if (debug) {
        malnBlk_dump(jb->joined, stderr, debugDumpOpts, "middle before common");
    }
    copyUnsharedGuideColumns(jb, jb->dests2, jb->cursor2, 0, jb->aln2CommonStart);
    if (debug) {
        malnBlk_dump(jb->joined, stderr, debugDumpOpts, "after before common");
    }

    // common
    if (jb->guideCommonStart < jb->guideCommonEnd) {
        assert(jb->dests1[0]->end == jb->guideCommonStart);
        joinSharedGuideColumns(jb);
        assert(jb->dests1[0]->end == jb->guideCommonEnd);
        if (debug) {
            malnBlk_dump(jb->joined, stderr, debugDumpOpts, "after common");
        }
    }

    // after common end
    copyUnsharedGuideColumns(jb, jb->dests1, jb->cursor1, jb->aln1CommonEnd, jb->blk1->alnWidth);
    if (debug) {
        malnBlk_dump(jb->joined, stderr, debugDumpOpts, "middle after common");
    }
    copyUnsharedGuideColumns(jb, jb->dests2, jb->cursor2, jb->aln2CommonEnd, jb->blk2->alnWidth);
    if (debug) {
        malnBlk_dump(jb->joined, stderr, debugDumpOpts, "after after common");
    }

    malnBlk_setTree(jb->joined, joinTrees(jb->cursor1, jb->cursor2, jb->srcDestCompMap));

    assertJoinedGuideComp(jb->dests1[0], jb->guide1, jb->guide2);
    assertJoinedComps(jb->cursor1, jb->dests1);
    assertJoinedComps(jb->cursor2, jb->dests2);
    malnBlk_finish(jb->joined);

    struct malnBlk *joined = jb->joined;
    malnJoinBlks_destruct(jb);
    return joined;
}
