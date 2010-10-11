#include "malnBlk.h"
#include "malnSet.h"
#include "malnComp.h"
#include "mafTree.h"
#include "genome.h"
#include "common.h"

/* constructor */
struct malnBlk *malnBlk_construct(void) {
    struct malnBlk *blk;
    AllocVar(blk);
    return blk;
}

/* set the tree on a block */
void malnBlk_setTree(struct malnBlk *blk, mafTree *mTree) {
    blk->mTree = mTree;
}

/* constructor clone */
struct malnBlk *malnBlk_constructClone(struct malnBlk *srcBlk) {
    struct malnBlk *blk = malnBlk_construct();
    for (struct malnComp *srcComp = srcBlk->comps; srcComp != NULL; srcComp = srcComp->next) {
        malnBlk_addComp(blk, malnComp_constructClone(srcComp));
    }
    slReverse(&blk->comps);
    malnBlk_setTree(blk, mafTree_clone(srcBlk->mTree, blk));
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
    assert(comp1->ncLink->treeOrder != comp2->ncLink->treeOrder);
    return comp1->ncLink->treeOrder - comp2->ncLink->treeOrder;
}

/* sort components by tree */
static void malnBlk_sortComps(struct malnBlk *blk) {
    assert(blk->mTree != NULL);
    mafTree_sortChildren(blk->mTree);
    cmpMTree = blk->mTree;
    slSort(&blk->comps, compTreeOrderCmp);
    cmpMTree = NULL;
}

/* finish construction a block, setting component attributes and sorting
 * components */
void malnBlk_finish(struct malnBlk *blk) {
    malnBlk_sortComps(blk); // FIXME: don't need to do this for all cases
    malnBlk_validate(blk); // produces better error messages
    malnBlk_assert(blk);   // really only to catch bugs in this code, not input
}

/* unlink a comp from the list. */
static void unlinkComp(struct malnBlk *blk, struct malnComp *comp) {
    struct malnComp *prevComp = NULL, *scanComp = blk->comps;
    while ((scanComp != comp) && (scanComp != NULL)) {
        prevComp = scanComp;
        scanComp = scanComp->next;
    }
    assert(scanComp != NULL);
    if (prevComp == NULL) {
        blk->comps = comp->next;
    } else {
        prevComp->next = comp->next;
    }
}

/* Unlink a component from the block, also removing it from the tree and
 * malnSet map, if its in a set.  Component is not freed. */
void malnBlk_unlink(struct malnBlk *blk, struct malnComp *comp) {
    unlinkComp(blk, comp);
    if (blk->malnSet != NULL) {
        malnSet_removeComp(blk->malnSet, comp);
    }
}

/* get the root component */
struct malnComp *malnBlk_getRootComp(struct malnBlk *blk) {
    struct malnComp *root = slLastEl(blk->comps);
    if ((malnComp_getLoc(root) & mafTreeLocRoot) == 0) {
        malnBlk_dump(blk, "last component is not root", stderr);
        errAbort("last component is not root");
    }
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
    struct malnBlk *rcBlk = malnBlk_construct();
    for (struct malnComp *srcComp = srcBlk->comps; srcComp != NULL; srcComp = srcComp->next) {
        malnBlk_addComp(rcBlk, malnComp_reverseComplement(srcComp));
    }
    slReverse(&rcBlk->comps);
    malnBlk_setTree(rcBlk, mafTree_clone(srcBlk->mTree, rcBlk));
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

/* check for consistency within a components of a MAF, generating an
 * error if there are not consistent */
void malnBlk_validate(struct malnBlk *blk) {
    // check for overlapping root components
    struct malnComp *rootComp = malnBlk_getRootComp(blk);
    for (struct malnComp *comp2 = blk->comps; comp2 != NULL; comp2 = comp2->next) {
        if ((comp2 != rootComp) && malnComp_overlap(rootComp, comp2)) {
            if (true) { // FIXME: tmp debug
                malnBlk_dump(blk, "overlapping components", stderr);
            }
            errAbort("overlapping root components detected with in a block: %s:%d-%d (%c) and %s:%d-%d (%c)",
                     rootComp->seq->orgSeqName, rootComp->start, rootComp->end, rootComp->strand,
                     comp2->seq->orgSeqName, comp2->start, comp2->end, comp2->strand);
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
#if 0 // FIXME, I think this is bogus for tandem dups
        // check for overlaps
        for (struct malnComp *comp2 = comp->next; comp2 != NULL; comp2 = comp2->next) {
            if (malnComp_overlap(comp, comp2)) { // FIXME: tmp
                malnBlk_dump(blk, "hash overlapping", stderr);
            }
            assert(!malnComp_overlap(comp, comp2));
        }
#endif
    }
#endif
}

/* compare two blocks for deterministic sorting */
int malnBlk_cmp(struct malnBlk *blk1, struct malnBlk *blk2) {
    int diff = malnComp_cmp(malnBlk_getRootComp(blk1), malnBlk_getRootComp(blk2));
    if (diff == 0) {
        diff = slCount(blk1->comps) - slCount(blk2->comps);
        if (diff == 0) {
            for (struct malnComp *comp1 = blk1->comps, *comp2 = blk2->comps; (diff == 0) && (comp1 != NULL); comp1 = comp1->next, comp2 = comp2->next) {
                diff = malnComp_cmp(comp1, comp2);
            }
        }
    }
    return diff;
}

/* print a block for debugging purposes */
void malnBlk_dump(struct malnBlk *blk, const char *label, FILE *fh) {
    char *nhTree = (blk->mTree != NULL) ? mafTree_format(blk->mTree) : cloneString("NULL");
    fprintf(fh, "%s %zx %d: %s\n", label, (size_t)blk, blk->alnWidth, nhTree);
    freeMem(nhTree);
    for (struct malnComp *comp = blk->comps; comp != NULL; comp = comp->next) {
        malnComp_dump(comp, "    ", fh);
    }    
}
