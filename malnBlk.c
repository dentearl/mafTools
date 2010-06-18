#include "malnBlk.h"
#include "malnComp.h"
#include "mafTree.h"
#include "common.h"

/* constructor */
struct malnBlk *malnBlk_construct(mafTree *mTree) {
    struct malnBlk *blk;
    AllocVar(blk);
    blk->mTree = mTree;
    return blk;
}

/* constructor clone */
struct malnBlk *malnBlk_constructClone(struct malnBlk *srcBlk) {
    struct malnBlk *blk = malnBlk_construct(srcBlk->mTree);
    for (struct malnComp *srcComp = srcBlk->comps; srcComp != NULL; srcComp = srcComp->next) {
        malnBlk_addComp(blk, malnComp_constructClone(srcComp));
    }
    malnBlk_sortComps(blk);
    return blk;
}

/* destructor */
void malnBlk_destruct(struct malnBlk *blk) {
    if (blk != NULL) {
        while (blk->comps != NULL) {
            malnComp_destruct(slPopHead(&blk->comps));
        }
        // FIXME: free tree?
        freeMem(blk);
    }
}

/* set location and type attributes for a component tree */
static void malnBlk_setCompLocAttr(struct malnBlk *blk, struct malnComp *comp, struct Genome *refGenome) {
    comp->isReference = (comp->seq->genome == refGenome);
    switch (mafTree_getLoc(blk->mTree, comp->seq->orgSeqName, comp->chromStart, comp->chromEnd)) {
    case mafTreeRoot:
        comp->treeLoc = malnCompTreeRoot;
        break;
    case mafTreeInternal:
        comp->treeLoc = malnCompTreeInternal;
        break;
    case mafTreeLeaf:
        comp->treeLoc = malnCompTreeLeaf;
        break;
    }
}

/* set location and type attributes from tree */
void malnBlk_setLocAttr(struct malnBlk *blk, struct Genome *refGenome) {
    for (struct malnComp *comp = blk->comps; comp != NULL; comp = comp->next) {
        malnBlk_setCompLocAttr(blk, comp, refGenome);
    }
}

/* add a component */
void malnBlk_addComp(struct malnBlk *blk, struct malnComp *comp) {
    comp->blk = blk;
    slAddHead(&blk->comps, comp);
    blk->alnWidth = max(blk->alnWidth, malnComp_getWidth(comp));
}

/* tree from block being compared */
static mafTree *cmpMTree = NULL;

/* compare function for two components objects by tree order */
static int compTreeOrderCmp(const void *vcomp1, const void *vcomp2) {
    const struct malnComp *comp1 = *((const struct malnComp**)vcomp1);
    const struct malnComp *comp2 = *((const struct malnComp**)vcomp2);
    return mafTree_treeOrderCmp(cmpMTree, comp1->seq->orgSeqName, comp1->chromStart, comp1->chromEnd, comp2->seq->orgSeqName, comp2->chromStart, comp2->chromEnd);
}

/* sort components by tree */
void malnBlk_sortComps(struct malnBlk *blk) {
    assert(blk->mTree != NULL);
    cmpMTree = blk->mTree;
    slSort(&blk->comps, compTreeOrderCmp);
    cmpMTree = NULL;
}

/* get the root component */
struct malnComp *malnBlk_getRootComp(struct malnBlk *blk) {
    struct malnComp *root = slLastEl(blk->comps);
    assert((root->treeLoc & malnCompTreeRoot) == 0);
    return root;
}

/* find a component by seq and start, NULL if not found  */
struct malnComp *malnBlk_findBySeqStart(struct malnBlk *blk, struct Seq *seq, int start) {
    for (struct malnComp *comp = blk->comps; comp != NULL; comp = comp->next) {
        if ((comp->seq == seq) && (comp->start == start)) {
            return comp;
        }
    }
    return NULL;
}

/* find a component by seq and chrom range, NULL if not found  */
struct malnComp *malnBlk_findCompByChromRange(struct malnBlk *blk, struct Seq *seq, int chromStart, int chromEnd) {
    for (struct malnComp *comp = blk->comps; comp != NULL; comp = comp->next) {
        if ((comp->seq == seq) && (comp->chromStart == chromStart) && (comp->chromEnd == chromEnd)) {
            return comp;
        }
    }
    return NULL;
}

/* block reverse complement */
struct malnBlk *malnBlk_reverseComplement(struct malnBlk *blk) {
    struct malnBlk *rcBlk = malnBlk_construct(blk->mTree);
    for (struct malnComp *comp = blk->comps; comp != NULL; comp = comp->next) {
        malnBlk_addComp(rcBlk, malnComp_reverseComplement(comp));
    }
    malnBlk_sortComps(rcBlk);
    return rcBlk;
}

/* Pad out an alignment so all columns are equal sized */
void malnBlk_pad(struct malnBlk *blk) {
    for (struct malnComp *comp = blk->comps; comp != NULL; comp = comp->next) {
        assert(malnComp_getWidth(comp) <= blk->alnWidth);
        if (malnComp_getWidth(comp) <= blk->alnWidth) {
            malnComp_pad(comp, blk->alnWidth);
        }
    }
}

/* assert that the block is set-consistent */
void malnBlk_assert(struct malnBlk *blk) {
#ifndef NDEBUG
   for (struct malnComp *comp = blk->comps; comp != NULL; comp = comp->next) {
        assert(malnComp_getWidth(comp) == blk->alnWidth);
   }
#endif
}
