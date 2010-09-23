#include "malnSet.h"
#include "malnComp.h"
#include "malnBlk.h"
#include "mafTree.h"
#include "sonLibSortedSet.h"
#include "sonLibList.h"
#include "common.h"
#include "jkmaf.h"
#include "genomeRangeTree.h"

struct malnSet {
    struct Genomes *genomes;
    stSortedSet *blks;
    struct genomeRangeTree *compRangeMap; /* range index slRefs of malnComp
                                           * objects.  chrom name is org.seq,
                                           * allowing indexing all components, not
                                           * just reference, need to join dups */
    bool dying;   /* indicates set is in the process of being deleted; used
                   * to prevent circular deletes */
};

/* Add all components to range map. */
static void addCompsToMap(struct malnSet *malnSet, struct malnBlk *blk) {
    for (struct malnComp *comp = blk->comps; comp != NULL; comp = comp->next) {
        genomeRangeTreeAddValList(malnSet->compRangeMap, comp->seq->orgSeqName, comp->chromStart, comp->chromEnd, slRefNew(comp));
    }
}

/* build the range tree when needed */
static void buildRangeTree(struct malnSet *malnSet) {
    malnSet->compRangeMap = genomeRangeTreeNew();
    stSortedSetIterator *iter = stSortedSet_getIterator(malnSet->blks);
    struct malnBlk *blk;
    while ((blk = stSortedSet_getNext(iter)) != NULL) {
        addCompsToMap(malnSet, blk);
    }
    stSortedSet_destructIterator(iter);
}

/* remove a single object from the range tree (missing genomeRangeTree functionality).
 * this just set reference to NULL, which must be handled when getting overlaps*/
static void removeCompRange(struct malnSet *malnSet, struct malnBlk *blk, struct malnComp *comp) {
    bool found = false;
    for (struct range *rng = genomeRangeTreeAllOverlapping(malnSet->compRangeMap, comp->seq->orgSeqName, comp->chromStart, comp->chromEnd); (rng != NULL) && (!found); rng = rng->next) {
        for (struct slRef *compRef = rng->val; (compRef != NULL) && (!found); compRef = compRef->next) {
            struct malnComp *comp = compRef->val;
            if ((comp !=NULL) && (comp->blk == blk)) {
                compRef->val = NULL;
                found = true;
            }
        }
    }
    assert(found);
}

/* remove all references to this blk from the rangeTree */
static void removeBlkRanges(struct malnSet *malnSet, struct malnBlk *blk) {
    for (struct malnComp *comp = blk->comps; comp != NULL; comp = comp->next) {
        removeCompRange(malnSet, blk, comp);
    }
}

/* convert a mafComp to an malnComp */
static struct malnComp *mafCompToMAlnComp(struct Genomes *genomes, struct mafComp *comp) {
    char buf[128];
    char *srcDb = mafCompGetSrcDb(comp, buf, sizeof(buf));
    if (srcDb == NULL) {
        errAbort("Error: no org name in MAF component, source must be org.seq: %s", comp->src);
    }
    return malnComp_construct(genomesObtainSeq(genomes, srcDb, mafCompGetSrcName(comp), comp->srcSize),
                              comp->strand, comp->start, comp->start+comp->size, comp->text);
}

/* convert a mafAli to an malnBlk */
static struct malnBlk *mafAliToMAlnBlk(struct Genomes *genomes, struct mafAli *ali, double defaultBranchLength, struct Genome *treelessRootGenome) {
    struct malnBlk *blk = malnBlk_construct(mafTree_constructFromMaf(genomes, ali, defaultBranchLength, treelessRootGenome));
    for (struct mafComp *comp = ali->components; comp != NULL; comp = comp->next) {
        malnBlk_addComp(blk, mafCompToMAlnComp(genomes, comp));
    }
    malnBlk_finish(blk);
#if 0
    malnBlk_assert(blk);
#endif
    return blk;
}

/* get associated genomes object  */
struct Genomes *malnSet_getGenomes(struct malnSet *malnSet) {
    return malnSet->genomes;
}

/* add a block to a malnSet */
void malnSet_addBlk(struct malnSet *malnSet, struct malnBlk *blk) {
    assert(blk->malnSet == NULL);
    malnBlk_assert(blk);
    blk->malnSet = malnSet;
    stSortedSet_insert(malnSet->blks, blk);
    if (malnSet->compRangeMap != NULL) {
        addCompsToMap(malnSet, blk);
    }
}

/* remove a block from malnSet */
void malnSet_removeBlk(struct malnSet *malnSet, struct malnBlk *blk) {
    assert(blk->malnSet == malnSet);
    if (!malnSet->dying) {
        if (malnSet->compRangeMap != NULL) {
            removeBlkRanges(malnSet, blk);
        }
        stSortedSet_remove(malnSet->blks, blk);
    }
    blk->malnSet = NULL;
}

/* remove a block from malnSet and free the block */
void malnSet_deleteBlk(struct malnSet *malnSet, struct malnBlk *blk) {
    malnSet_removeBlk(malnSet, blk);
    malnBlk_destruct(blk);
}

/* construct an empty malnSet  */
struct malnSet *malnSet_construct(struct Genomes *genomes) {
    struct malnSet *malnSet;
    AllocVar(malnSet);
    malnSet->genomes = genomes;
    malnSet->blks = stSortedSet_construct2((void (*)(void *))malnBlk_destruct);
    return malnSet;
}

/* Construct a malnSet from a MAF file. defaultBranchLength is used to
 * assign branch lengths when inferring trees from the MAF. */
struct malnSet *malnSet_constructFromMaf(struct Genomes *genomes, char *mafFileName, double defaultBranchLength, struct Genome *treelessRootGenome) {
    struct malnSet *malnSet = malnSet_construct(genomes);
    struct mafFile *mafFile = mafOpen(mafFileName);
    struct mafAli *ali;
    while ((ali = mafNext(mafFile)) != NULL) {
        malnSet_addBlk(malnSet, mafAliToMAlnBlk(genomes, ali, defaultBranchLength, treelessRootGenome));
        mafAliFree(&ali);
    }
    malnSet_assert(malnSet);
    mafFileFree(&mafFile);
    return malnSet;
}

/* free slRef objects in the compRangeMap */
static void destructCompRangeMap(struct malnSet *malnSet) {
    struct hashCookie cookie = hashFirst(malnSet->compRangeMap->hash);
    struct hashEl *hel;
    while ((hel = hashNext(&cookie)) != NULL) {
        struct rbTree *rangeTree = hel->val;
        for (struct range *rng = rangeTreeList(rangeTree); rng != NULL; rng = rng->next) {
            slFreeList(&rng->val);
        }
    }
    genomeRangeTreeFree(&malnSet->compRangeMap);
}

/* destructor */
void malnSet_destruct(struct malnSet *malnSet) {
    // prevent attempts to remove from stSortedSet while stSortedSet
    // is being deleted.
    malnSet->dying = true;
    if (malnSet->compRangeMap != NULL) {
        destructCompRangeMap(malnSet);
    }
    stSortedSet_destruct(malnSet->blks);
    freeMem(malnSet);
}

/* convert a malnComp to a mafComp */
static struct mafComp *malnCompToMafComp(struct malnComp *comp) {
    struct mafComp *mc;
    AllocVar(mc);
    mc->src = cloneString(comp->seq->orgSeqName);
    mc->srcSize = comp->seq->size;
    mc->strand = comp->strand;
    mc->start = comp->start;
    mc->size = comp->end - comp->start;
    mc->text = cloneString(malnComp_getAln(comp));
    return mc;
}

/* convert a malnBlk to a mafAli */
static struct mafAli *malnAliToMafAli(struct malnBlk *blk) {
    malnBlk_assert(blk);
    struct mafAli *ma;
    AllocVar(ma);
    if (blk->mTree != NULL) {
        ma->tree = mafTree_format(blk->mTree);
    }
    ma->textSize = blk->alnWidth;
    for (struct malnComp *comp = blk->comps; comp != NULL; comp = comp->next) {
        slAddHead(&ma->components, malnCompToMafComp(comp));
    }
    slReverse(&ma->components);
    return ma;
}

/* compare the root components on two blocks */
static int blkCmpRootComp(const void *vblk1, const void *vblk2) {
    struct malnBlk *blk1 = (struct malnBlk *)vblk1;
    struct malnBlk *blk2 = (struct malnBlk *)vblk2;
    struct malnComp *root1 = slLastEl(blk1->comps);
    struct malnComp *root2 = slLastEl(blk2->comps);
    int diff = strcmp(root1->seq->genome->name, root2->seq->genome->name);
    if (diff == 0) {
        diff = strcmp(root1->seq->name, root2->seq->name);
    }
    if (diff == 0) {
        diff = root1->chromStart - root2->chromStart;
    }
    if (diff == 0) {
        diff = root1->chromEnd - root2->chromEnd;
    }
    return diff;
}

/* build a list of blocks, sorted by the root components */
static stList *buildRootSorted(struct malnSet *malnSet) {
    stList *sorted = stList_construct();
    stSortedSetIterator *iter = stSortedSet_getIterator(malnSet->blks);
    struct malnBlk *blk;
    while ((blk = stSortedSet_getNext(iter)) != NULL) {
        stList_append(sorted, blk);
    }
    stSortedSet_destructIterator(iter);
    stList_sort(sorted, blkCmpRootComp);
    return sorted;
}

/* write a block to a MAF */
static void writeBlkToMaf(struct malnBlk *blk, FILE *mafFh) {
    malnBlk_validate(blk);
    struct mafAli *ma = malnAliToMafAli(blk);
    mafWrite(mafFh, ma);
    mafAliFree(&ma);
}

/* write a malnSet to a MAF file  */
void malnSet_writeMaf(struct malnSet *malnSet, char *mafFileName) {
    malnSet_assert(malnSet);
    stList *sorted = buildRootSorted(malnSet);

    FILE *mafFh = mustOpen(mafFileName, "w");
    mafWriteStart(mafFh, NULL);
    stListIterator *iter = stList_getIterator(sorted);
    struct malnBlk *blk;
    while ((blk = stList_getNext(iter)) != NULL) {
        writeBlkToMaf(blk, mafFh);
    }
    stList_destructIterator(iter);
    mafWriteEnd(mafFh);
    carefulClose(&mafFh);
    stList_destruct(sorted);
}

/* get iterator of the blocks. Don't remove or add blocks while in motion. */
stSortedSetIterator *malnSet_getBlocks(struct malnSet *malnSet) {
    return stSortedSet_getIterator(malnSet->blks);
}

/* Get a list of components that overlaps the specified reference range
 * and are in blocks not flagged as done and passing treeLoc filters.
 * Return NULL if no overlaps. */
stList *malnSet_getOverlappingPendingComps(struct malnSet *malnSet, struct Seq *seq, int chromStart, int chromEnd, unsigned treeLocFilter) {
    if (malnSet->compRangeMap == NULL) {
        buildRangeTree(malnSet);
    }
    stList *overBlks = NULL;
    for (struct range *rng = genomeRangeTreeAllOverlapping(malnSet->compRangeMap, seq->orgSeqName, chromStart, chromEnd); rng != NULL; rng = rng->next) {
        for (struct slRef *compRef = rng->val; compRef != NULL; compRef = compRef->next) {
            struct malnComp *comp = compRef->val;
            // FIXME: not sure why the overlap check is needed, shouldn't return non-overlaping,
            // but it did!
            if ((comp != NULL) && ((comp->treeLoc & treeLocFilter) != 0) && (!comp->blk->done)
                && malnComp_overlapRange(comp, seq, chromStart, chromEnd)) {
                if (overBlks == NULL) { 
                    overBlks = stList_construct();
                }
                stList_append(overBlks, comp);
            }
        }
    }
    return overBlks;
}

/* assert some sanity checks on a set */
void malnSet_assert(struct malnSet *malnSet) {
#ifndef NDEBUG
    stSortedSetIterator *iter = stSortedSet_getIterator(malnSet->blks);
    struct malnBlk *blk;
    while ((blk = stSortedSet_getNext(iter)) != NULL) {
        assert(blk->malnSet == malnSet);
        malnBlk_assert(blk);
    }
    stSortedSet_destructIterator(iter);
#endif
}

/* assert done flag is set on all blocks */
void malnSet_assertDone(struct malnSet *malnSet) {
#ifndef NDEBUG
    stSortedSetIterator *iter = stSortedSet_getIterator(malnSet->blks);
    struct malnBlk *blk;
    while ((blk = stSortedSet_getNext(iter)) != NULL) {
        assert(blk->done);
    }
    stSortedSet_destructIterator(iter);
#endif
}

/* clear done flag on all blocks */
void malnSet_clearDone(struct malnSet *malnSet) {
    stSortedSetIterator *iter = stSortedSet_getIterator(malnSet->blks);
    struct malnBlk *blk;
    while ((blk = stSortedSet_getNext(iter)) != NULL) {
        blk->done = false;
    }
    stSortedSet_destructIterator(iter);
}

