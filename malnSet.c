#include "malnSet.h"
#include "malnComp.h"
#include "malnBlk.h"
#include "mafTree.h"
#include "malnBlkSet.h"
#include "sonLibList.h"
#include "sonLibString.h"
#include "stSafeC.h"
#include "common.h"
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
    malnBlk_assert(blk);
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
        malnBlk_assert(blk);
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
        malnBlk_dump(blk, fh, "  blk");
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
