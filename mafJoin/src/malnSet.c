#include "malnSet.h"
#include "malnComp.h"
#include "malnBlk.h"
#include "mafTree.h"
#include "malnBlkSet.h"
#include "malnCompCursor.h"
#include "sonLibList.h"
#include "sonLibString.h"
#include "stSafeC.h"
#include "common.h"
#include "jkmaf.h"
#include "genomeRangeTree.h"

struct malnSet {
    struct Genomes *genomes;
    char *mafFileName;       // maybe NULL
    struct malnBlkSet *blks;
    struct genomeRangeTree *compRangeMap; /* range index slRefs of malnComp
                                           * objects.  chrom name is org.seq,
                                           * allowing indexing all components, not
                                           * just guide, need to join dups */
    bool dying;   /* indicates set is in the process of being deleted; used
                   * to prevent circular deletes */
    struct malnBlkSet *dyingBlks;  // bocks flaged for deletion
};

/* add a component to the range map. */
void malnSet_addComp(struct malnSet *malnSet, struct malnComp *comp) {
    if (malnSet->compRangeMap != NULL) {
        genomeRangeTreeAddValList(malnSet->compRangeMap, comp->seq->orgSeqName, comp->chromStart, comp->chromEnd, slRefNew(comp));
    }
}

/* Add all components to range map. */
static void addCompsToMap(struct malnSet *malnSet, struct malnBlk *blk) {
    for (struct malnComp *comp = blk->comps; comp != NULL; comp = comp->next) {
        malnSet_addComp(malnSet, comp);
    }
}

/* build the range tree when needed */
static void buildRangeTree(struct malnSet *malnSet) {
    malnSet->compRangeMap = genomeRangeTreeNew();
    struct malnBlkSetIterator *iter = malnBlkSet_getIterator(malnSet->blks);
    struct malnBlk *blk;
    while ((blk = malnBlkSetIterator_getNext(iter)) != NULL) {
        addCompsToMap(malnSet, blk);
    }
    malnBlkSetIterator_destruct(iter);
}

/* convert a mafComp to an malnComp */
static struct malnComp *mafCompToMAlnComp(struct Genomes *genomes, struct mafComp *mafComp) {
    if ((mafComp->srcSize < 0) || (mafComp->start < 0) || (mafComp->size < 0)
        || ((mafComp->start+mafComp->size) > mafComp->srcSize)
        || (!((mafComp->strand == '+') || (mafComp->strand == '-')))) {
        errAbort("invalid mafComp: %s srcSize=%d strand=%c, start=%d size=%d",
                 mafComp->src, mafComp->srcSize, mafComp->strand, mafComp->start, mafComp->size);
    }
        
    char buf[128];
    char *srcDb = mafCompGetSrcDb(mafComp, buf, sizeof(buf));
    if (srcDb == NULL) {
        errAbort("Error: no org name in MAF component, source must be org.seq: %s", mafComp->src);
    }
    struct malnComp *comp = malnComp_constructFromStr(genomesObtainSeq(genomes, srcDb, mafCompGetSrcName(mafComp), mafComp->srcSize),
                                                      mafComp->strand, mafComp->start, mafComp->start+mafComp->size, mafComp->text);
    
    return comp;
}

static struct Genome *sortTreelessRootGenome = NULL;  // used by treeless sort

/* compare components to sort the root genome components last, with the
 * longest component the root. */
static int orderTreelessCmp(const void *vcomp1, const void *vcomp2) {
    struct malnComp *comp1 = *((struct malnComp **)vcomp1);
    struct malnComp *comp2 = *((struct malnComp **)vcomp2);
    if ((comp1->seq->genome != sortTreelessRootGenome) && (comp2->seq->genome != sortTreelessRootGenome)) {
        return 0;
    }
    if ((comp1->seq->genome == sortTreelessRootGenome) && (comp2->seq->genome != sortTreelessRootGenome)) {
        return 1;
    }
    if ((comp1->seq->genome != sortTreelessRootGenome) && (comp2->seq->genome == sortTreelessRootGenome)) {
        return -1;
    }
    return -(malnComp_getNumAligned(comp1) - malnComp_getNumAligned(comp2));
}

/* set the ordering for block were tree must be constructed */
static void orderTreeless(struct malnBlk *blk, struct Genome *treelessRootGenome) {
    sortTreelessRootGenome = treelessRootGenome;
    slSort(&(blk->comps), orderTreelessCmp);
    sortTreelessRootGenome = NULL;
}

/* convert a mafAli to an malnBlk */
static struct malnBlk *mafAliToMAlnBlk(struct Genomes *genomes, struct mafAli *mafAli, double defaultBranchLength, struct Genome *treelessRootGenome) {
    struct malnBlk *blk = malnBlk_construct();
    for (struct mafComp *mafComp = mafAli->components; mafComp != NULL; mafComp = mafComp->next) {
        malnBlk_addComp(blk, mafCompToMAlnComp(genomes, mafComp));
    }
    if (mafAli->tree != NULL) {
        slReverse(&blk->comps);
        malnBlk_setTree(blk, mafTree_constructFromNewick(genomes, blk, mafAli->tree));
    } else if (treelessRootGenome != NULL) {
        orderTreeless(blk, treelessRootGenome);
        malnBlk_setTree(blk, mafTree_constructFromAlignment(genomes, blk, defaultBranchLength));
    } else {
        errAbort("no tree in mafAli block and no treelessRoot genome specified");
    }
    malnBlk_finish(blk);
    return blk;
}

/* get associated genomes object  */
struct Genomes *malnSet_getGenomes(struct malnSet *malnSet) {
    return malnSet->genomes;
}

/* get the source MAF file or NULL if not specified */
char *malnSet_getMafFileName(struct malnSet *malnSet) {
    return malnSet->mafFileName;
}

/* add a block to a malnSet */
void malnSet_addBlk(struct malnSet *malnSet, struct malnBlk *blk) {
    assert(blk->malnSet == NULL);
    malnBlk_assert(blk, true);
    blk->malnSet = malnSet;
    malnBlkSet_add(malnSet->blks, blk);
    if (malnSet->compRangeMap != NULL) {
        addCompsToMap(malnSet, blk);
    }
}
 
/* add all blocks in to a malnSet */
void malnSet_addBlks(struct malnSet *malnSet, struct malnBlkSet *blks) {
    struct malnBlkSetIterator *iter = malnBlkSet_getIterator(blks);
    struct malnBlk *blk;
    while ((blk = malnBlkSetIterator_getNext(iter)) != NULL) {
        malnSet_addBlk(malnSet, blk);
    }
    malnBlkSetIterator_destruct(iter);
}


/* remove a single object from the range tree (missing genomeRangeTree functionality).
 * this just set guide to NULL, which must be handled when getting overlaps*/
static void removeCompRange(struct malnSet *malnSet, struct malnComp *comp) {
    bool found = false;
    for (struct range *rng = genomeRangeTreeAllOverlapping(malnSet->compRangeMap, comp->seq->orgSeqName, comp->chromStart, comp->chromEnd); (rng != NULL) && (!found); rng = rng->next) {
        for (struct slRef *compRef = rng->val; (compRef != NULL) && (!found); compRef = compRef->next) {
            if (compRef->val == comp) {
                compRef->val = NULL;
                found = true;
            }
        }
    }
    assert(found);
}

/* remove a single component from malnSet range map */
void malnSet_removeComp(struct malnSet *malnSet, struct malnComp *comp) {
    if (malnSet->compRangeMap != NULL) {
        removeCompRange(malnSet, comp);
    }
}

/* remove all guides to this blk from the rangeTree */
static void removeBlkRanges(struct malnSet *malnSet, struct malnBlk *blk) {
    for (struct malnComp *comp = blk->comps; comp != NULL; comp = comp->next) {
        removeCompRange(malnSet, comp);
    }
    malnBlkSet_remove(malnSet->blks, blk);
    blk->malnSet = NULL;
}

/* remove a block from malnSet */
void malnSet_removeBlk(struct malnSet *malnSet, struct malnBlk *blk) {
    assert(blk->malnSet == malnSet);
    if (!malnSet->dying) {
        if (malnSet->compRangeMap != NULL) {
            removeBlkRanges(malnSet, blk);
        }
        malnBlkSet_remove(malnSet->blks, blk);
    }
    blk->malnSet = NULL;
}

/* pop a block from the set.  order block are returned is arbitrary */
struct malnBlk *malnSet_popBlk(struct malnSet *malnSet) {
    struct malnBlk *blk = malnBlkSet_pop(malnSet->blks);
    if (blk != NULL) {
        if ((malnSet->compRangeMap != NULL) && !malnSet->dying) {
            removeBlkRanges(malnSet, blk);
        }
        blk->malnSet = NULL;
    }
    return blk;
}

/* remove a block from malnSet and free the block */
static void malnSet_deleteBlk(struct malnSet *malnSet, struct malnBlk *blk) {
    assert(blk->malnSet == malnSet);
    malnSet_removeBlk(malnSet, blk);
    malnBlk_destruct(blk);
}

/* construct an empty malnSet.  mafFileName maybe NULL  */
struct malnSet *malnSet_construct(struct Genomes *genomes, char *mafFileName) {
    struct malnSet *malnSet;
    AllocVar(malnSet);
    malnSet->genomes = genomes;
    if (mafFileName != NULL) {
        malnSet->mafFileName = cloneString(mafFileName);
    }
    // n.b. don't try to sort on a key, just use addess key, otherwise
    // deleting entries in merge will delete the wrong entries.
    malnSet->blks = malnBlkSet_construct();
    malnSet->dyingBlks = malnBlkSet_construct();
    return malnSet;
}

/* Initialize the alignment starts in an array of splits. Adjust bounds of
 * partitions so that the have least some bases of the block's root aligned,
 * as require by split code.  Return actual number of partitions */
static int splitBlockInitStarts(struct malnBlk *blk, int numParts, int partSize, int alnStarts[]) {
    struct malnComp *rootComp = malnBlk_getRootComp(blk);
    int alnNext = 0;
    int iPart = 0;
    // fill in combining boundary as needed to make sure there are
    // some root bases.
    while ((iPart < numParts) && (alnNext < blk->alnWidth)) {
        int alnEnd = min(alnNext+partSize, blk->alnWidth);
        if (malnComp_anyAlignedRange(rootComp, alnNext, alnEnd)) {
            alnStarts[iPart++] = alnNext;
        }
        alnNext = alnEnd;
    }
    alnStarts[iPart] = alnNext;

    // if the last entry does have any aligned root bases, move the bound back
    while ((iPart > 0) && !malnComp_anyAlignedRange(rootComp, alnStarts[iPart-1], alnStarts[iPart])) {
        iPart--;
        alnStarts[iPart-1] = alnStarts[iPart];
    }
    return iPart;
}

/* split a block in to roughly equal sizes pieces under the limit and add
 * them.  This is very annoying as each piece must include at least some
 * bases of the root sequence. */
static void splitBlock(struct malnSet *malnSet, int maxInputBlkWidth, struct malnBlk *blk) {
    int numParts = (blk->alnWidth + maxInputBlkWidth - 1) / maxInputBlkWidth;
    int partSize = (blk->alnWidth + numParts - 1) / numParts;
    int alnStarts[numParts+1];
    int usedNumParts = splitBlockInitStarts(blk, numParts, partSize, alnStarts);

    for (int iPart = 0; iPart < usedNumParts; iPart++) {
        malnSet_addBlk(malnSet, malnBlk_constructSubrange(blk, alnStarts[iPart], alnStarts[iPart+1]));
    }
}

/* add a mafAli to the set */
static void addMafAli(struct malnSet *malnSet, struct mafAli *mafAli, int maxInputBlkWidth, double defaultBranchLength, struct Genome *treelessRootGenome) {
    struct malnBlk *blk = mafAliToMAlnBlk(malnSet->genomes, mafAli, defaultBranchLength, treelessRootGenome);
    if (blk->alnWidth < maxInputBlkWidth) {
        malnSet_addBlk(malnSet, blk);
    } else {
        splitBlock(malnSet, maxInputBlkWidth, blk);
        malnBlk_destruct(blk);
    }
}

/* Construct a malnSet from a MAF file. defaultBranchLength is used to
 * assign branch lengths when inferring trees from the MAF. */
struct malnSet *malnSet_constructFromMaf(struct Genomes *genomes, char *mafFileName, int maxInputBlkWidth, double defaultBranchLength, struct Genome *treelessRootGenome) {
    struct malnSet *malnSet = malnSet_construct(genomes, mafFileName);
    struct mafFile *mafFile = mafOpen(mafFileName);
    struct mafAli *ali;
    while ((ali = mafNext(mafFile)) != NULL) {
        addMafAli(malnSet, ali, maxInputBlkWidth, defaultBranchLength, treelessRootGenome);
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

/* delete all blocks in a set */
static void deleteAllBlks(struct malnSet *malnSet, struct malnBlkSet *blkMap) {
    struct malnBlk *blk;
    while ((blk = malnBlkSet_pop(blkMap)) != NULL) {
        if (blk->malnSet == NULL) {
            malnBlk_destruct(blk);
        } else {
            malnSet_deleteBlk(malnSet, blk);
        }
    }
}

/* destructor */
void malnSet_destruct(struct malnSet *malnSet) {
    // prevent attempts to remove from stSortedSet while stSortedSet
    // is being deleted.
    malnSet->dying = true;
    if (malnSet->compRangeMap != NULL) {
        destructCompRangeMap(malnSet);
    }
    deleteAllBlks(malnSet, malnSet->blks);
    malnBlkSet_destruct(malnSet->blks);
    malnBlkSet_destruct(malnSet->dyingBlks);
    freeMem(malnSet->mafFileName);
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
    // FIXME: do block at a time
    mc->text = needMem(comp->alnWidth+1);
    struct malnCompCursor cursor = malnCompCursor_make(comp);
    while (malnCompCursor_incr(&cursor)) {
        mc->text[cursor.alnIdx] = malnCompCursor_getCol(&cursor);
    }
    return mc;
}

/* convert a malnBlk to a mafAli */
static struct mafAli *malnAliToMafAli(struct malnBlk *blk) {
    malnBlk_assert(blk, true);
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
    struct malnComp *root1 = malnBlk_getRootComp(blk1);
    struct malnComp *root2 = malnBlk_getRootComp(blk2);
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
    struct malnBlkSetIterator *iter = malnBlkSet_getIterator(malnSet->blks);
    struct malnBlk *blk;
    while ((blk = malnBlkSetIterator_getNext(iter)) != NULL) {
        stList_append(sorted, blk);
    }
    malnBlkSetIterator_destruct(iter);
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
struct malnBlkSetIterator *malnSet_getBlocks(struct malnSet *malnSet) {
    return malnBlkSet_getIterator(malnSet->blks);
}

/* Get a table with the set of blocks. Useful when one needs to add remove or
 * add blocks while scanning. */
struct malnBlkSet *malnSet_getBlockSetCopy(struct malnSet *malnSet) {
    return malnBlkSet_constructClone(malnSet->blks);
}

/* compare function for component list */
static int sortCompListCmpFn(const void *a, const void *b) {
    struct malnComp *ca = (struct malnComp*)a;
    struct malnComp *cb = (struct malnComp*)b;
    int diff = ca->blk->alnWidth - cb->blk->alnWidth;
    if (diff == 0) {
        diff = malnComp_cmp(ca, cb); 
    }
    return diff;
}

/* check is a component in the overlap list should be included */
static bool keepOverlap(struct malnComp *comp, struct Seq *seq, int chromStart, int chromEnd, unsigned treeLocFilter) {
    // FIXME: not sure why the overlap check is needed, shouldn't return non-overlaping,
    // but it did!  Need to verify this.
    return ((comp != NULL) && ((malnComp_getLoc(comp) & treeLocFilter) != 0) && (!comp->blk->deleted) && malnComp_overlapRange(comp, seq, chromStart, chromEnd));
}        

/* Get a list of components that overlap the specified guide range and are in
 * blocks not flagged as dying and matches treeLoc filters.  Return NULL if no
 * overlaps.  List is sorted by ascending width, which helps the merge
 * efficiency.  */
stList *malnSet_getOverlappingComps(struct malnSet *malnSet, struct Seq *seq, int chromStart, int chromEnd, unsigned treeLocFilter) {
    if (malnSet->compRangeMap == NULL) {
        buildRangeTree(malnSet);
    }
    stList *overComps = NULL;
    for (struct range *rng = genomeRangeTreeAllOverlapping(malnSet->compRangeMap, seq->orgSeqName, chromStart, chromEnd); rng != NULL; rng = rng->next) {
        for (struct slRef *compRef = rng->val; compRef != NULL; compRef = compRef->next) {
            struct malnComp *comp = compRef->val;
            if (keepOverlap(comp, seq, chromStart, chromEnd, treeLocFilter)) {
                if (overComps == NULL) { 
                    overComps = stList_construct();
                }
                stList_append(overComps, comp);
            }
        }
    }

    // sort so tests are reproducible
    if (overComps != NULL) {
        stList_sort(overComps, sortCompListCmpFn);
    }
    return overComps;
}

/* Get a list of components that overlap or are adjacent to the specified
 * guide range and are in blocks not flagged as dying and matches treeLoc
 * filters.  Return NULL if no overlaps.  List is sorted by ascending width,
 * which helps the merge efficiency. */
stList *malnSet_getOverlappingAdjacentComps(struct malnSet *malnSet, struct Seq *seq, int chromStart, int chromEnd, unsigned treeLocFilter) {
    return malnSet_getOverlappingComps(malnSet, seq, chromStart-1, chromEnd+1, treeLocFilter);
}

/* assert some sanity checks on a set */
void malnSet_assert(struct malnSet *malnSet) {
#ifndef NDEBUG
    struct malnBlkSetIterator *iter = malnBlkSet_getIterator(malnSet->blks);
    struct malnBlk *blk;
    while ((blk = malnBlkSetIterator_getNext(iter)) != NULL) {
        assert(blk->malnSet == malnSet);
        malnBlk_assert(blk, true);
    }
    malnBlkSetIterator_destruct(iter);
#endif
}

/* Record a block as deleted.  It is allowed to add blocks with are not
 * associated with a malnSet, as a way of cleaning up intermediate blocks.
 */
void malnSet_markAsDeleted(struct malnSet *malnSet, struct malnBlk *blk) {
    assert((blk->malnSet == malnSet) || (blk->malnSet == NULL));
    blk->deleted = true;
    malnBlkSet_add(malnSet->dyingBlks, blk);
    malnBlk_freeSeqMem(blk);
}

/* get the fraction of blocks that are dying */
float malnSet_fractionDying(struct malnSet *malnSet) {
    float total = malnBlkSet_size(malnSet->blks);
    if (total == 0.0) {
        return 0.0;
    } else {
        return ((float)malnBlkSet_size(malnSet->dyingBlks)) / total;
    }
}

/* delete blocks marked as dying */
void malnSet_deleteDying(struct malnSet *malnSet) {
    deleteAllBlks(malnSet, malnSet->dyingBlks);
}

/* print set for debugging */
void malnSet_dumpv(struct malnSet *malnSet, FILE *fh, const char *label, va_list args) {
    char *fmtLabel = stSafeCDynFmtv(label, args);
    fprintf(fh, "malnSet: begin %s\n", fmtLabel);
    struct malnBlkSetIterator *iter = malnBlkSet_getIterator(malnSet->blks);
    struct malnBlk *blk;
    while ((blk = malnBlkSetIterator_getNext(iter)) != NULL) {
        fputs("  ", fh);
        malnBlk_dump(blk, fh, malnBlk_dumpDefault, "  blk");
    }
    malnBlkSetIterator_destruct(iter);
    fprintf(fh, "malnSet: end %s\n", fmtLabel);
    freeMem(fmtLabel);
}

/* print set for debugging */
void malnSet_dump(struct malnSet *malnSet, FILE *fh, const char *label, ...) {
    va_list args;
    va_start(args, label);
    malnSet_dumpv(malnSet, fh, label, args);
    va_end(args);
}

/* print set for debugging to a file */
void malnSet_dumpFile(struct malnSet *malnSet, char *dumpFile, const char *label, ...) {
    FILE *fh = mustOpen(dumpFile, "w");
    va_list args;
    va_start(args, label);
    malnSet_dumpv(malnSet, fh, label, args);
    va_end(args);
    carefulClose(&fh);
}

/* debug dump to a file in a directory.  Will do nothing if dumpDir is NULL */
void malnSet_dumpToDir(struct malnSet *malnSet, char *dumpDir, char *setName, char *stepName) {
    if (dumpDir != NULL) {
        char label[PATH_LEN], dumpFile[PATH_LEN];
        safef(label, sizeof(label), "%s-%s", setName, stepName);
        safef(dumpFile, sizeof(dumpFile), "%s/%s-%s.dmp", dumpDir, setName, stepName);
        malnSet_dumpFile(malnSet, dumpFile, "%s", label);
    }
}
