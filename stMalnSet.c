#include "stMalnSet.h"
#include "stMalnComp.h"
#include "stMalnBlk.h"
#include "stMafTree.h"
#include "common.h"
#include "jkmaf.h"
#include "genomeRangeTree.h"

struct stMalnSet {
    struct Genomes *genomes;
    struct Genome *refGenome;
    struct stMalnBlk *blks;
    struct genomeRangeTree *refBlks; /* range index of reference genome components to slRefs to stMalnBlk
                                      * objects.  Built on demand. */
};

/* get associated genomes object  */
struct Genomes *stMalnSet_getGenomes(struct stMalnSet *malnSet) {
    return malnSet->genomes;
}

/* convert a mafComp to an stMalnComp */
static struct stMalnComp *mafCompToMAlnComp(struct Genomes *genomes, struct mafComp *comp) {
    char buf[128];
    char *srcDb = mafCompGetSrcDb(comp, buf, sizeof(buf));
    if (srcDb == NULL) {
        errAbort("Error: no org name in MAF component, source must be org.seq: %s", comp->src);
    }
    return stMalnComp_construct(genomesObtainSeq(genomes, srcDb, mafCompGetSrcName(comp), comp->srcSize),
                                comp->strand, comp->start, comp->start+comp->size, comp->text);
}

/* convert a mafAli to an stMalnBlk */
static struct stMalnBlk *mafAliToMAlnBlk(struct Genomes *genomes, struct mafAli *ali, double defaultBranchLength) {
    struct stMalnBlk *blk = stMalnBlk_construct(stMafTree_constructFromMaf(ali, defaultBranchLength));
    for (struct mafComp *comp = ali->components; comp != NULL; comp = comp->next) {
        stMalnBlk_addComp(blk, mafCompToMAlnComp(genomes, comp));
    }
    stMalnBlk_sortComps(blk);
    return blk;
}

/* add a block to a malnSet */
void stMalnSet_addBlk(struct stMalnSet *malnSet, struct stMalnBlk *blk) {
    stMalnBlk_assert(blk);
    slAddHead(&malnSet->blks, blk);
}

/* construct an empty stMalnSet  */
struct stMalnSet *stMalnSet_construct(struct Genomes *genomes, struct Genome *refGenome) {
    struct stMalnSet *malnSet;
    AllocVar(malnSet);
    malnSet->genomes = genomes;
    malnSet->refGenome = refGenome;
    return malnSet;
}

/* Construct a stMalnSet from a MAF file. defaultBranchLength is used to
 * assign branch lengths when inferring trees from pair-wise MAFs. */
struct stMalnSet *stMalnSet_constructFromMaf(struct Genomes *genomes, struct Genome *refGenome, char *mafFileName, double defaultBranchLength) {
    struct stMalnSet *malnSet = stMalnSet_construct(genomes, refGenome);
    struct mafFile *mafFile = mafOpen(mafFileName);
    struct mafAli *ali;
    while ((ali = mafNext(mafFile)) != NULL) {
        stMalnSet_addBlk(malnSet, mafAliToMAlnBlk(genomes, ali, defaultBranchLength));
        mafAliFree(&ali);
    }
    slReverse(&malnSet->blks);
    return malnSet;
}

/* convert a stMalnComp to a mafComp */
static struct mafComp *malnCompToMafComp(struct stMalnComp *comp) {
    struct mafComp *mc;
    AllocVar(mc);
    mc->src = cloneString(comp->seq->orgSeqName);
    mc->srcSize = comp->seq->size;
    mc->strand = comp->strand;
    mc->start = comp->start;
    mc->size = comp->end - comp->start;
    mc->text = cloneString(stMalnComp_getAln(comp));
    return mc;
}

/* convert a stMalnBlk to a mafAli */
static struct mafAli *malnAliToMafAli(struct stMalnBlk *blk) {
    stMalnBlk_assert(blk);
    struct mafAli *ma;
    AllocVar(ma);
    if (blk->mTree != NULL) {
        ma->tree = stMafTree_format(blk->mTree);
    }
    ma->textSize = blk->alnWidth;
    for (struct stMalnComp *comp = blk->comps; comp != NULL; comp = comp->next) {
        slAddHead(&ma->components, malnCompToMafComp(comp));
    }
    slReverse(&ma->components);
    return ma;
}

/* write a block to a MAF */
static void writeBlkToMaf(struct stMalnBlk *blk, FILE *mafFh) {
    struct mafAli *ma = malnAliToMafAli(blk);
    mafWrite(mafFh, ma);
    mafAliFree(&ma);
}

/* write a stMalnSet to a MAF file  */
void stMalnSet_writeMaf(struct stMalnSet *malnSet, char *mafFileName) {
    FILE *mafFh = mustOpen(mafFileName, "w");
    mafWriteStart(mafFh, NULL);
    for (struct stMalnBlk *blk = malnSet->blks; blk != NULL; blk = blk->next) {
        writeBlkToMaf(blk, mafFh);
    }
    mafWriteEnd(mafFh);
    carefulClose(&mafFh);
}

/* get list of blocks.  Don't modify this list. */
struct stMalnBlk *stMalnSet_getBlocks(struct stMalnSet *malnSet) {
    return malnSet->blks;
}

/* add block to range map for all reference genome components */
static void addRefRanges(struct stMalnSet *malnSet, struct stMalnBlk *blk) {
    for (struct stMalnComp *mc = blk->comps; mc != NULL; mc = mc->next) {
        if (mc->seq->genome == malnSet->refGenome) {
            genomeRangeTreeAddValList(malnSet->refBlks, mc->seq->name, mc->chromStart, mc->chromEnd, slRefNew(blk));
        }
    }
}

/* build the range tree when needed */
static void buildRangeTree(struct stMalnSet *malnSet) {
    malnSet->refBlks = genomeRangeTreeNew();
    for (struct stMalnBlk *blk = malnSet->blks; blk != NULL; blk = blk->next) {
        addRefRanges(malnSet, blk);
    }
}

/* get list of slRefs to blks who's reference range overlaps the specified
 * range.  slFreeList of results when done. */
struct slRef *stMalnSet_getOverlapping(struct stMalnSet *malnSet, struct Seq *seq, int chromStart, int chromEnd) {
    if (malnSet->refBlks == NULL) {
        buildRangeTree(malnSet);
    }
    struct slRef *head = NULL;
    for (struct range *rng = genomeRangeTreeAllOverlapping(malnSet->refBlks, seq->name, chromStart, chromEnd); rng != NULL; rng = rng->next) {
        struct slRef *blkRef = rng->val;
        slAddHead(&head, slRefNew(blkRef->val));
    }
    return head;
}
