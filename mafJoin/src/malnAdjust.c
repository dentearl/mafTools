#include "common.h"
#include "malnSet.h"
#include "malnBlkSet.h"
#include "malnBlk.h"


/* make width adjustment to a block. */
static void adjustBlkWidth(struct malnBlk *inBlk, struct malnSet *malnSetOut, int maxBlkWidth) {
    int alnStart = 0;
    while (alnStart < inBlk->alnWidth) {
        int alnEnd = alnStart + maxBlkWidth;
        if (alnEnd > inBlk->alnWidth) {
            alnEnd = inBlk->alnWidth;
        }
        malnSet_addBlk(malnSetOut, malnBlk_constructSubrange(inBlk, alnStart, alnEnd));
        alnStart = alnEnd;
    }
}

/* make adjustments to blocks in a set, adding to a new set */
void malnAdjust_adjustSet(struct malnSet *malnSetIn, struct malnSet *malnSetOut, int maxBlkWidth) {
    struct malnBlkSetIterator *iter = malnSet_getBlocks(malnSetIn);
    struct malnBlk *inBlk;
    while ((inBlk = malnBlkSetIterator_getNext(iter)) != NULL) {
        adjustBlkWidth(inBlk, malnSetOut, maxBlkWidth);
    }
    malnBlkSetIterator_destruct(iter);
}

