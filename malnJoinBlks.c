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
    struct malnBlkCursor *cursor1;   // cursors into input block
    struct malnBlkCursor *cursor2;
    struct malnBlk *joined;          // joined block
    struct malnComp **dests1;        // arrays of source to destination components, first common
    struct malnComp **dests2;
    struct malnCompCompMap *srcDestCompMap; // map of source to destination components for tree join
    struct malnBlk *freeBlk;   // block to free if not NULL (due to reverse-complement);
};

/* create new reference components from the two being joined */
static struct malnComp *createJoinedRefComp(struct malnJoinBlks *jb) {
    int start = min(jb->ref1->start, jb->ref2->start);
    struct malnComp *comp = malnComp_construct(jb->ref1->seq, jb->ref1->strand, start, start, "");
    malnBlk_addComp(jb->joined, comp);
    malnCompCompMap_insert(jb->srcDestCompMap, jb->ref1, comp);
    malnCompCompMap_insert(jb->srcDestCompMap, jb->ref2, comp);
    return comp;
}

/* add non-reference components to the joined blk and destination component
 * array */
static void addCompsToJoined(struct malnJoinBlks *jb, struct malnBlkCursor *blkCursor, struct malnComp **destComps) {
    // skip rows[0], which is the reference
    for (int i = 1; i < blkCursor->numRows; i++) {
        destComps[i] = malnComp_construct(blkCursor->rows[i].comp->seq, blkCursor->rows[i].comp->strand, blkCursor->rows[i].comp->start, blkCursor->rows[i].comp->start, "");
        malnBlk_addComp(jb->joined, destComps[i]);
        malnCompCompMap_insert(jb->srcDestCompMap, blkCursor->rows[i].comp, destComps[i]);
    }
}

/* construct malnJoinBlks state object for the join */
static struct malnJoinBlks *malnJoinBlks_construct(struct malnComp *refComp1, struct malnComp *refComp2) {
    struct malnJoinBlks *jb;
    AllocVar(jb);
    jb->ref1 = refComp1;
    jb->ref2 = refComp2;
    jb->blk1 = refComp1->blk;
    jb->blk2 = refComp2->blk;
    if (debug) { // FIXME: tmp
        malnComp_dump(jb->ref1, "refComp1", stderr);
        malnBlk_dump(jb->blk1, "inBlk1", stderr);
        malnComp_dump(jb->ref2, "refComp2", stderr);
        malnBlk_dump(jb->blk2, "inBlk2", stderr);
    }

    // reverse complement if needed (must do before making cursors)
    if (jb->ref1->strand != jb->ref2->strand) {
        jb->blk2 = jb->freeBlk = malnBlk_reverseComplement(jb->blk2);
        jb->ref2 = malnBlk_findCompByChromRange(jb->blk2, jb->ref2->seq, jb->ref2->chromStart, jb->ref2->chromEnd);
        assert(jb->ref2 != NULL);
        if (debug) { // FIXME: tmp
            malnBlk_dump(jb->blk2, "inBlk2rc", stderr);
        }
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
    jb->blk1->done = true;  // must do before destruct, as blk1 might be freeBlk
    jb->blk2->done = true;
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
static void copyUnsharedRefColumns(struct malnBlk *blkJoined, struct malnComp **destComps, struct malnBlkCursor *blkCursor, int alnStart, int alnEnd) {
    malnBlk_assert(blkJoined);  // FIXME
    malnBlk_assert(blkCursor->blk);  // FIXME
    if (debug) {
        malnBlk_dump(blkJoined, "blkJoined", stderr);  // FIXME
        malnBlk_dump(blkCursor->blk, "blkCursor", stderr);  // FIXME
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
    struct malnJoinBlks *jb = malnJoinBlks_construct(refComp1, refComp2);

    int refCommonStart = max(jb->ref1->start, jb->ref2->start);
    int refCommonEnd = min(jb->ref1->end, jb->ref2->end);
    int aln1CommonStart, aln1CommonEnd, aln2CommonStart, aln2CommonEnd;
    if (!(malnComp_seqRangeToAlnRange(jb->ref1, refCommonStart, refCommonEnd, &aln1CommonStart, &aln1CommonEnd)
          && malnComp_seqRangeToAlnRange(jb->ref2, refCommonStart, refCommonEnd, &aln2CommonStart, &aln2CommonEnd))) {
        errAbort("BUG: failure to get alignment ranges for common reference sequence range");
    }
    if (joinedCompRet != NULL) {
        *joinedCompRet = jb->dests1[0];
    }

    if (debug) {
        fprintf(stderr, "refCommon: %s:%d-%d\n", jb->ref1->seq->orgSeqName, refCommonStart, refCommonEnd);
        fprintf(stderr, "align1Commom:: %d-%d\n", aln1CommonStart, aln1CommonEnd);
        fprintf(stderr, "align2Commom:: %d-%d\n", aln2CommonStart, aln2CommonEnd);
    }

    // before common start
    copyUnsharedRefColumns(jb->joined, jb->dests1, jb->cursor1, 0, aln1CommonStart);
    if (debug) {malnBlk_dump(jb->joined, "joined@1", stderr);}
    copyUnsharedRefColumns(jb->joined, jb->dests2, jb->cursor2, 0, aln2CommonStart);
    if (debug) {malnBlk_dump(jb->joined, "joined@2", stderr);}

    // common
    assert(jb->dests1[0]->end == refCommonStart);
    joinSharedRefColumns(jb->joined, jb->dests1, jb->cursor1, aln1CommonEnd, jb->dests2, jb->cursor2, aln2CommonEnd);
    assert(jb->dests1[0]->end == refCommonEnd);
    if (debug) {malnBlk_dump(jb->joined, "joined@3", stderr);}

    // after common end
    copyUnsharedRefColumns(jb->joined, jb->dests1, jb->cursor1, aln1CommonEnd, jb->blk1->alnWidth);
    copyUnsharedRefColumns(jb->joined, jb->dests2, jb->cursor2, aln2CommonEnd, jb->blk2->alnWidth);
    
    malnBlk_setTree(jb->joined, joinTrees(jb->cursor1, jb->cursor2, jb->srcDestCompMap));

    assertJoinedRefComp(jb->dests1[0], jb->ref1, jb->ref2);
    assertJoinedComps(jb->cursor1, jb->dests1);
    assertJoinedComps(jb->cursor2, jb->dests2);
    malnBlk_finish(jb->joined);
    malnBlk_assert(jb->joined);

    struct malnBlk *joined = jb->joined;
    malnJoinBlks_destruct(jb);
    return joined;
}
