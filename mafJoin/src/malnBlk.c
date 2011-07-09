#include "malnBlk.h"
#include "malnSet.h"
#include "malnComp.h"
#include "malnCompCursor.h"
#include "mafTree.h"
#include "malnCompCompMap.h"
#include "stSafeC.h"
#include "genome.h"
#include "common.h"

/* constructor */
struct malnBlk *malnBlk_construct(void) {
    static int nextObjId = 1;
    struct malnBlk *blk;
    AllocVar(blk);
    blk->objId = nextObjId++;
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
    blk->alnWidth = srcBlk->alnWidth;
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

/* If a block is part of a set, mark it for deletion in that set,
 * otherwise deleted immediately*/
void malnBlk_markOrDelete(struct malnBlk *blk) {
    if (blk->malnSet != NULL) {
        malnSet_markAsDeleted(blk->malnSet, blk);
    } else {
        malnBlk_destruct(blk);
    }
}

/* Clear the sequences for a dying block.  This supports freeing up large
 * amounts of memory immediately without having to update the metadata,
 * requiring iterators to be invalidated, etc. */
void malnBlk_freeSeqMem(struct malnBlk *blk) {
    // FIXME: might be better just to set blk/comp pointers to null in metadata,
    // but this is easy for now
    assert(blk->deleted);
    for (struct malnComp *comp = blk->comps; comp != NULL; comp = comp->next) {
        malnComp_freeSeqMem(comp);
    }
}

/* add a component */
void malnBlk_addComp(struct malnBlk *blk, struct malnComp *comp) {
    comp->blk = blk;
    slAddHead(&blk->comps, comp);
    blk->alnWidth = max(blk->alnWidth, comp->alnWidth);
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
    malnBlk_validate(blk);  // produces better error messages
#ifdef ASSERT_SLOW
    malnBlk_assert(blk, false);    // really only to catch bugs in this code, not input
#endif
}

/* Unlink a component from the block, also removing it from the tree and
 * malnSet map, if its in a set.  Component is not freed. */
void malnBlk_unlink(struct malnBlk *blk, struct malnComp *comp) {
    assert(comp->blk == blk);
    slRemoveEl(&blk->comps, comp);
    mafTree_deleteNode(blk->mTree, comp->ncLink);
    if (blk->malnSet != NULL) {
        malnSet_removeComp(blk->malnSet, comp);
    }
    comp->blk = NULL;
}

/* get the root component */
struct malnComp *malnBlk_getRootComp(struct malnBlk *blk) {
    return slLastEl(blk->comps);
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
    malnBlk_assert(srcBlk, true);
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
        assert(comp->alnWidth <= blk->alnWidth);
        if (comp->alnWidth <= blk->alnWidth) {
            malnComp_pad(comp, blk->alnWidth);
        }
    }
}

/* check for consistency within a components of a MAF, generating an
 * error if there are not consistent */
void malnBlk_validate(struct malnBlk *blk) {
    // check sizes
    for (struct malnComp *comp = blk->comps; comp != NULL; comp = comp->next) {
        if (comp->alnWidth != blk->alnWidth) {
            errAbort("component width (%d) doesn't match alignment width (%d): %s:%d-%d",
                     comp->alnWidth, blk->alnWidth,
                     comp->seq->orgSeqName, comp->chromStart, comp->chromEnd);
        }
    }

    // check for overlapping root components
    struct malnComp *rootComp = malnBlk_getRootComp(blk);
    for (struct malnComp *comp2 = blk->comps; comp2 != NULL; comp2 = comp2->next) {
        if ((comp2 != rootComp) && malnComp_overlap(rootComp, comp2)) {
            errAbort("overlapping root components detected with in a block: %s:%d-%d (%c) and %s:%d-%d (%c)",
                     rootComp->seq->orgSeqName, rootComp->start, rootComp->end, rootComp->strand,
                     comp2->seq->orgSeqName, comp2->start, comp2->end, comp2->strand);
        }
    }
}

#ifdef ASSERT_SLOW
/* assert that there are no columns that have nothing aligned in them */
static void assertNoEmptyColumns(struct malnBlk *blk) {
    bool colsWithBases[blk->alnWidth];
    memset(colsWithBases, false, blk->alnWidth);
    for (struct malnComp *comp = blk->comps; comp != NULL; comp = comp->next) {
        for (struct malnCompSeg *seg = comp->segs; seg != NULL; seg = seg->next) {
            memset(colsWithBases+seg->alnStart, true, seg->size);
        }
    }
    bool *notSet = memchr(colsWithBases, false, blk->alnWidth);
    if (notSet != NULL) {
        malnBlk_dump(blk, stderr, malnBlk_dumpDefault, "empty column at %ld", (notSet - colsWithBases));
    }
    assert(notSet == NULL);
}
#endif

/* assert that the block is set-consistent.  ncLinkCheck can be false for
 * sanity check before links and tree are constructed. */
void malnBlk_assert(struct malnBlk *blk, bool ncLinkCheck) {
#ifndef NDEBUG
    assert(blk->comps != NULL);
    assert(blk->alnWidth >= 0);
    for (struct malnComp *comp = blk->comps; comp != NULL; comp = comp->next) {
        assert(comp->blk == blk);
        if (comp->alnWidth != blk->alnWidth) {
            malnBlk_dump(blk, stderr, malnBlk_dumpDefault, "invalid component width");
        }
        assert(comp->alnWidth == blk->alnWidth);
        malnComp_assert(comp, ncLinkCheck);
    }
    assertNoEmptyColumns(blk);
    if (ncLinkCheck) {
        mafTree_assert(blk->mTree, blk);
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

/* compare blocks by ascending width */
int malnBlk_assendingWidthCmp(struct malnBlk *blk1, struct malnBlk *blk2) {
    return blk1->alnWidth - blk2->alnWidth;
}

/* take subrange of a component */
static void compSubRange(struct malnBlk *subBlk, struct malnComp *comp, int alnStart, int alnEnd, struct malnCompCompMap *srcDestCompMap) {
    struct malnComp *subComp = malnComp_constructSubrange(comp, alnStart, alnEnd);
    if (subComp != NULL) {
        malnBlk_addComp(subBlk, subComp);
    }
    malnCompCompMap_add(srcDestCompMap, comp, subComp);
}

/* construct an alignment block from a subrange of this block */
struct malnBlk *malnBlk_constructSubrange(struct malnBlk *blk, int alnStart, int alnEnd) {
    assert((alnStart < alnEnd) && (alnStart >= 0) && (alnEnd <= blk->alnWidth));
    struct malnBlk *subBlk = malnBlk_construct();
    struct malnCompCompMap *srcDestCompMap = malnCompCompMap_construct();
    for (struct malnComp *comp = blk->comps; comp != NULL; comp = comp->next) {
        compSubRange(subBlk, comp, alnStart, alnEnd, srcDestCompMap);
    }
    subBlk->mTree = mafTree_subrangeClone(blk->mTree, srcDestCompMap);
    malnBlk_finish(subBlk);
    malnCompCompMap_destruct(srcDestCompMap);
    return subBlk;
}

/* print a block for debugging purposes */
void malnBlk_dumpv(struct malnBlk *blk, FILE *fh, unsigned opts, const char *label, va_list args) {
    char *nhTree = (blk->mTree != NULL) ? mafTree_format(blk->mTree) : cloneString("NULL");
    char *fmtLabel = stSafeCDynFmtv(label, args);
    fprintf(fh, "%s #%d width: %d%s %s\n", fmtLabel, blk->objId, blk->alnWidth, (blk->deleted ? " deleted" : ""), nhTree);
    if ((opts & malnBlk_dumpInclTree) && (blk->mTree != NULL)) {
        mafTree_dump(blk->mTree, fh);
    }
    freeMem(fmtLabel);
    freeMem(nhTree);
    for (struct malnComp *comp = blk->comps; comp != NULL; comp = comp->next) {
        if (opts & malnBlk_dumpNoSeqs) {
            fputc('\t', fh);
            malnComp_prInfo(comp, fh);
            fputc('\n', fh);
        } else {
            malnComp_dump(comp, fh, "\t");
        }
    }    
}
/* print a block for debugging purposes */
void malnBlk_dump(struct malnBlk *blk, FILE *fh, unsigned opts, const char *label, ...) {
    va_list args;
    va_start(args, label);
    malnBlk_dumpv(blk, fh, opts, label, args);
    va_end(args);
}
