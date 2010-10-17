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

static bool debug = false;  // FIXME: tmp

/* Object use to keep state */
struct malnJoinBlks {
    struct malnComp *ref1;           // reference components
    struct malnComp *ref2;
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
    int refCommonStart;  // common reference sequence coordinates, zero length for adjacent
    int refCommonEnd;
    int aln1CommonStart; // mapping of the ref coordinates to input alignments
    int aln1CommonEnd;
    int aln2CommonStart;
    int aln2CommonEnd;
};

/* compute alignment coordinates */
static void calcAlignmentCoords(struct malnJoinBlks *jb) {
    jb->refCommonStart = max(jb->ref1->start, jb->ref2->start);
    jb->refCommonEnd = min(jb->ref1->end, jb->ref2->end);
    if (jb->refCommonStart < jb->refCommonEnd) {
        // overlapping
        if (!(malnComp_seqRangeToAlnRange(jb->ref1, jb->refCommonStart, jb->refCommonEnd, &jb->aln1CommonStart, &jb->aln1CommonEnd)
              && malnComp_seqRangeToAlnRange(jb->ref2, jb->refCommonStart, jb->refCommonEnd, &jb->aln2CommonStart, &jb->aln2CommonEnd))) {
            errAbort("BUG: failure to get alignment ranges for common reference sequence range");
        }
    } else {
        // adjacent
        jb->aln1CommonStart = jb->aln1CommonEnd = jb->cursor1->alnWidth;
        jb->aln2CommonStart = jb->aln2CommonEnd = jb->cursor2->alnWidth;
    }
}

/* create new reference components from the two being joined */
static struct malnComp *createJoinedRefComp(struct malnJoinBlks *jb) {
    int start = min(jb->ref1->start, jb->ref2->start);
    struct malnComp *comp = malnComp_construct(jb->ref1->seq, jb->ref1->strand, start, start, "");
    malnBlk_addComp(jb->joined, comp);
    malnCompCompMap_add(jb->srcDestCompMap, jb->ref1, comp);
    malnCompCompMap_add(jb->srcDestCompMap, jb->ref2, comp);
    return comp;
}

/* add non-reference components to the joined blk and destination component
 * array */
static void addCompsToJoined(struct malnJoinBlks *jb, struct malnBlkCursor *blkCursor, struct malnComp **destComps) {
    // skip rows[0], which is the reference
    for (int i = 1; i < blkCursor->numRows; i++) {
        destComps[i] = malnComp_construct(blkCursor->rows[i].comp->seq, blkCursor->rows[i].comp->strand, blkCursor->rows[i].comp->start, blkCursor->rows[i].comp->start, "");
        malnBlk_addComp(jb->joined, destComps[i]);
        malnCompCompMap_add(jb->srcDestCompMap, blkCursor->rows[i].comp, destComps[i]);
    }
}

/* construct malnJoinBlks state object for the join */
static struct malnJoinBlks *malnJoinBlks_construct(struct malnComp *refComp1, struct malnComp *refComp2) {
    struct malnJoinBlks *jb;
    AllocVar(jb);
    jb->ref1 = refComp1;
    jb->ref2 = refComp2;
    jb->blk1 = refComp1->blk;
    jb->blk2 = jb->inBlk2 = refComp2->blk;

    // reverse complement if needed (must do before making cursors)
    if (jb->ref1->strand != jb->ref2->strand) {
        jb->blk2 = jb->freeBlk = malnBlk_reverseComplement(jb->blk2);
        jb->ref2 = malnBlk_findCompByChromRange(jb->blk2, jb->ref2->seq, jb->ref2->chromStart, jb->ref2->chromEnd);
        assert(jb->ref2 != NULL);
    }

    jb->cursor1 = malnBlkCursor_construct(jb->blk1, jb->ref1, NULL);
    jb->cursor2 = malnBlkCursor_construct(jb->blk2, jb->ref2, NULL);
    jb->joined = malnBlk_construct();
    jb->dests1 = needMem(jb->cursor1->numRows * sizeof(struct malnComp *));
    jb->dests2 = needMem(jb->cursor2->numRows * sizeof(struct malnComp *));
    jb->srcDestCompMap = malnCompCompMap_construct();
    jb->freeBlk = NULL;

    jb->dests1[0] = jb->dests2[0] = createJoinedRefComp(jb);
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
static void copyUnsharedRefColumns(struct malnJoinBlks *jb, struct malnComp **destComps, struct malnBlkCursor *blkCursor, int alnStart, int alnEnd) {
  for (int i = 0; i < blkCursor->numRows; i++) {
        malnComp_appendCompAln(destComps[i], blkCursor->rows[i].comp, alnStart, alnEnd);
    }
    malnBlkCursor_setAlignCol(blkCursor, alnEnd);
    jb->joined->alnWidth += (alnEnd - alnStart);  // FIXME should have append methods
    malnBlk_pad(jb->joined);
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
static void copySharedRefColumns(struct malnJoinBlks *jb) {
    assert(isRefAligned(jb->cursor1));
    assert(isRefAligned(jb->cursor2));
    while (isRefAligned(jb->cursor1) && isRefAligned(jb->cursor2) && (jb->cursor1->alnIdx < jb->aln1CommonEnd) && (jb->cursor2->alnIdx < jb->aln2CommonEnd)) {
        copyColumn(jb->dests1, jb->cursor1, false);
        copyColumn(jb->dests2, jb->cursor2, true);
        jb->joined->alnWidth++;  // FIXME should have append methods
    }
}

/* copy contiguous unaligned-to-reference columns to join blk */
static void copyUnalignedSharedColumns(struct malnBlk *blkJoined, struct malnComp **destComps, struct malnBlkCursor *blkCursor, int alnCommonEnd) {
    while ((!isRefAligned(blkCursor)) && (blkCursor->alnIdx < alnCommonEnd)) {
        copyColumn(destComps, blkCursor, false);
        blkJoined->alnWidth++;  // FIXME should have append methods
    }
    malnBlk_pad(blkJoined);
}

/* join columns based on shared reference sequence regions */
static void joinSharedRefColumns(struct malnJoinBlks *jb) {
    assert(jb->cursor1->rows[0].pos == jb->cursor2->rows[0].pos);
    while ((jb->cursor1->alnIdx < jb->aln1CommonEnd) && (jb->cursor1->alnIdx < jb->aln1CommonEnd)) {
        copySharedRefColumns(jb);
        copyUnalignedSharedColumns(jb->joined, jb->dests1, jb->cursor1, jb->aln1CommonEnd);
        copyUnalignedSharedColumns(jb->joined, jb->dests2, jb->cursor2, jb->aln2CommonEnd);
    }
    assert(jb->cursor1->alnIdx == jb->aln1CommonEnd);
    assert(jb->cursor2->alnIdx == jb->aln2CommonEnd);
}

/* join two blocks using their specified reference components.  Optionally return
 * resulting join component. */
struct malnBlk *malnJoinBlks(struct malnComp *refComp1, struct malnComp *refComp2, struct malnComp **joinedCompRet) {
    assert(malnComp_overlapAdjacent(refComp1, refComp2));
    struct malnJoinBlks *jb = malnJoinBlks_construct(refComp1, refComp2);
    if (joinedCompRet != NULL) {
        *joinedCompRet = jb->dests1[0];
    }
    calcAlignmentCoords(jb);
    if (debug) { // FIXME: tmp
        fprintf(stderr, "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n");
        malnComp_dump(jb->ref1, "refComp1", stderr);
        malnBlk_dump(jb->blk1, "inBlk1", stderr);
        malnComp_dump(jb->ref2, "refComp2", stderr);
        malnBlk_dump(jb->inBlk2, "inBlk2", stderr);
        if (jb->blk2 != jb->inBlk2) {
            malnBlk_dump(jb->blk2, "inBlk2rc", stderr);
        }
    }

    // before common start
    copyUnsharedRefColumns(jb, jb->dests1, jb->cursor1, 0, jb->aln1CommonStart);
    copyUnsharedRefColumns(jb, jb->dests2, jb->cursor2, 0, jb->aln2CommonStart);

    // common
    if (jb->refCommonStart < jb->refCommonEnd) {
        assert(jb->dests1[0]->end == jb->refCommonStart);
        joinSharedRefColumns(jb);
        assert(jb->dests1[0]->end == jb->refCommonEnd);
    }

    // after common end
    copyUnsharedRefColumns(jb, jb->dests1, jb->cursor1, jb->aln1CommonEnd, jb->blk1->alnWidth);
    copyUnsharedRefColumns(jb, jb->dests2, jb->cursor2, jb->aln2CommonEnd, jb->blk2->alnWidth);
    
    malnBlk_setTree(jb->joined, joinTrees(jb->cursor1, jb->cursor2, jb->srcDestCompMap));

    assertJoinedRefComp(jb->dests1[0], jb->ref1, jb->ref2);
    assertJoinedComps(jb->cursor1, jb->dests1);
    assertJoinedComps(jb->cursor2, jb->dests2);
    malnBlk_finish(jb->joined);
    malnBlk_assert(jb->joined);

    struct malnBlk *joined = jb->joined;
    malnJoinBlks_destruct(jb);
    if (debug) { // FIXME: tmp
        malnBlk_dump(joined, "joined", stderr);
        fprintf(stderr, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
    }
    return joined;
}
