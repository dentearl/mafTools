#include "malnBlk.h"
#include "malnSet.h"
#include "malnComp.h"
#include "mafTree.h"
#include "genome.h"
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
    struct malnBlk *blk = malnBlk_construct(mafTree_clone(srcBlk->mTree));
    for (struct malnComp *srcComp = srcBlk->comps; srcComp != NULL; srcComp = srcComp->next) {
        malnBlk_addComp(blk, malnComp_constructClone(srcComp));
    }
    malnBlk_finish(blk);
    return blk;
}

/* destructor */
void malnBlk_destruct(struct malnBlk *blk) {
    if (blk != NULL) {
        if (blk->malnSet != NULL) {
            malnSet_removeBlk(blk->malnSet, blk);
        }
        while (blk->comps != NULL) {
            malnComp_destruct(slPopHead(&blk->comps));
        }
        mafTree_destruct(blk->mTree);
        freeMem(blk);
    }
}

/* set the tree location attribute for a component tree.  
 * FIXME: this feels hacky compared to saving links with tree and
 * check on the fly. */
static void malnBlk_setCompLocAttr(struct malnBlk *blk, struct malnComp *comp) {
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
static void malnBlk_setLocAttr(struct malnBlk *blk) {
    for (struct malnComp *comp = blk->comps; comp != NULL; comp = comp->next) {
        malnBlk_setCompLocAttr(blk, comp);
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
static void malnBlk_sortComps(struct malnBlk *blk) {
    assert(blk->mTree != NULL);
    cmpMTree = blk->mTree;
    slSort(&blk->comps, compTreeOrderCmp);
    cmpMTree = NULL;
}

/* finish construction a block, setting component attributes and sorting
 * components */
void malnBlk_finish(struct malnBlk *blk) {
    malnBlk_sortComps(blk);
    malnBlk_setLocAttr(blk);
    malnBlk_assert(blk);
}

/* get the root component */
struct malnComp *malnBlk_getRootComp(struct malnBlk *blk) {
    struct malnComp *root = slLastEl(blk->comps);
    assert((root->treeLoc & malnCompTreeRoot) != 0);
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
struct malnBlk *malnBlk_reverseComplement(struct malnBlk *srcBlk) {
    malnBlk_assert(srcBlk);
    struct malnBlk *rcBlk = malnBlk_construct(mafTree_clone(srcBlk->mTree));
    for (struct malnComp *srcComp = srcBlk->comps; srcComp != NULL; srcComp = srcComp->next) {
        malnBlk_addComp(rcBlk, malnComp_reverseComplement(srcComp));
    }
    malnBlk_finish(rcBlk);
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
    assert(blk->comps != NULL);
    for (struct malnComp *comp = blk->comps; comp != NULL; comp = comp->next) {
        assert(comp->blk == blk);
        assert(malnComp_getWidth(comp) == blk->alnWidth);
        malnComp_assert(comp);
    }
#endif
}

/* print a block for debugging purposes */
void malnBlk_dump(struct malnBlk *blk, const char *label, FILE *fh) {
    fprintf(fh, "%s %zx %d\n", label, (size_t)blk, blk->alnWidth);
    for (struct malnComp *comp = blk->comps; comp != NULL; comp = comp->next) {
        malnComp_dump(comp, "    ", fh);
    }    
}
