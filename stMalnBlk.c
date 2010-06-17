#include "stMalnBlk.h"
#include "stMalnComp.h"
#include "stMafTree.h"
#include "common.h"

/* constructor */
struct stMalnBlk *stMalnBlk_construct(stMafTree *mTree) {
    struct stMalnBlk *blk;
    AllocVar(blk);
    blk->mTree = mTree;
    return blk;
}

/* constructor clone */
struct stMalnBlk *stMalnBlk_constructClone(struct stMalnBlk *srcBlk) {
    struct stMalnBlk *blk = stMalnBlk_construct(srcBlk->mTree);
    for (struct stMalnComp *srcComp = srcBlk->comps; srcComp != NULL; srcComp = srcComp->next) {
        stMalnBlk_addComp(blk, stMalnComp_constructClone(srcComp));
    }
    stMalnBlk_sortComps(blk);
    return blk;
}

/* destructor */
void stMalnBlk_destruct(struct stMalnBlk *blk) {
    if (blk != NULL) {
        while (blk->comps != NULL) {
            stMalnComp_destruct(slPopHead(&blk->comps));
        }
        // FIXME: free tree?
        freeMem(blk);
    }
}

/* set location and type attributes for a component tree */
static void stMalnBlk_setCompLocAttr(struct stMalnBlk *blk, struct stMalnComp *comp, struct Genome *refGenome) {
    comp->isReference = (comp->seq->genome == refGenome);
    switch (stMafTree_getLoc(blk->mTree, comp->seq->orgSeqName, comp->chromStart, comp->chromEnd)) {
    case stMafTreeRoot:
        comp->treeLoc = stMalnCompTreeRoot;
        break;
    case stMafTreeInternal:
        comp->treeLoc = stMalnCompTreeInternal;
        break;
    case stMafTreeLeaf:
        comp->treeLoc = stMalnCompTreeLeaf;
        break;
    }
}

/* set location and type attributes from tree */
void stMalnBlk_setLocAttr(struct stMalnBlk *blk, struct Genome *refGenome) {
    for (struct stMalnComp *comp = blk->comps; comp != NULL; comp = comp->next) {
        stMalnBlk_setCompLocAttr(blk, comp, refGenome);
    }
}

/* add a component */
void stMalnBlk_addComp(struct stMalnBlk *blk, struct stMalnComp *comp) {
    slAddHead(&blk->comps, comp);
    blk->alnWidth = max(blk->alnWidth, stMalnComp_getWidth(comp));
}

/* tree from block being compared */
static stMafTree *cmpMTree = NULL;

/* compare function for two components objects by tree order */
static int compTreeOrderCmp(const void *vcomp1, const void *vcomp2) {
    const struct stMalnComp *comp1 = *((const struct stMalnComp**)vcomp1);
    const struct stMalnComp *comp2 = *((const struct stMalnComp**)vcomp2);
    return stMafTree_treeOrderCmp(cmpMTree, comp1->seq->orgSeqName, comp1->chromStart, comp1->chromEnd, comp2->seq->orgSeqName, comp2->chromStart, comp2->chromEnd);
}

/* sort components by tree */
void stMalnBlk_sortComps(struct stMalnBlk *blk) {
    assert(blk->mTree != NULL);
    cmpMTree = blk->mTree;
    slSort(&blk->comps, compTreeOrderCmp);
    cmpMTree = NULL;
}

/* find a component by seq and start, NULL if not found  */
struct stMalnComp *stMalnBlk_findBySeqStart(struct stMalnBlk *blk, struct Seq *seq, int start) {
    for (struct stMalnComp *comp = blk->comps; comp != NULL; comp = comp->next) {
        if ((comp->seq == seq) && (comp->start == start)) {
            return comp;
        }
    }
    return NULL;
}

/* find a component by seq and chrom range, NULL if not found  */
struct stMalnComp *stMalnBlk_findCompByChromRange(struct stMalnBlk *blk, struct Seq *seq, int chromStart, int chromEnd) {
    for (struct stMalnComp *comp = blk->comps; comp != NULL; comp = comp->next) {
        if ((comp->seq == seq) && (comp->chromStart == chromStart) && (comp->chromEnd == chromEnd)) {
            return comp;
        }
    }
    return NULL;
}

/* block reverse complement */
struct stMalnBlk *stMalnBlk_reverseComplement(struct stMalnBlk *blk) {
    struct stMalnBlk *rcBlk = stMalnBlk_construct(blk->mTree);
    for (struct stMalnComp *comp = blk->comps; comp != NULL; comp = comp->next) {
        stMalnBlk_addComp(rcBlk, stMalnComp_reverseComplement(comp));
    }
    stMalnBlk_sortComps(rcBlk);
    return rcBlk;
}

/* Pad out an alignment so all columns are equal sized */
void stMalnBlk_pad(struct stMalnBlk *blk) {
    for (struct stMalnComp *comp = blk->comps; comp != NULL; comp = comp->next) {
        assert(stMalnComp_getWidth(comp) <= blk->alnWidth);
        if (stMalnComp_getWidth(comp) <= blk->alnWidth) {
            stMalnComp_pad(comp, blk->alnWidth);
        }
    }
}

/* assert that the block is set-consistent */
void stMalnBlk_assert(struct stMalnBlk *blk) {
#ifndef NDEBUG
   for (struct stMalnComp *comp = blk->comps; comp != NULL; comp = comp->next) {
        assert(stMalnComp_getWidth(comp) == blk->alnWidth);
   }
#endif
}
