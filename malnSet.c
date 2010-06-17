#include "malnSet.h"
#include "malnComp.h"
#include "malnBlk.h"
#include "mafTree.h"
#include "common.h"
#include "jkmaf.h"
#include "genomeRangeTree.h"

struct malnSet {
    struct Genomes *genomes;
    struct Genome *refGenome;
    struct malnBlk *blks;
    struct genomeRangeTree *refBlks; /* range index of reference genome components to slRefs to malnBlk
                                      * objects.  Built on demand. */
};

/* get associated genomes object  */
struct Genomes *malnSet_getGenomes(struct malnSet *malnSet) {
    return malnSet->genomes;
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
static struct malnBlk *mafAliToMAlnBlk(struct Genomes *genomes, struct Genome *refGenome, struct mafAli *ali, double defaultBranchLength) {
    struct malnBlk *blk = malnBlk_construct(mafTree_constructFromMaf(ali, defaultBranchLength));
    for (struct mafComp *comp = ali->components; comp != NULL; comp = comp->next) {
        malnBlk_addComp(blk, mafCompToMAlnComp(genomes, comp));
    }
    malnBlk_sortComps(blk);
    malnBlk_setLocAttr(blk, refGenome);
    return blk;
}

/* add a block to a malnSet */
void malnSet_addBlk(struct malnSet *malnSet, struct malnBlk *blk) {
    malnBlk_assert(blk);
    slAddHead(&malnSet->blks, blk);
}

/* construct an empty malnSet  */
struct malnSet *malnSet_construct(struct Genomes *genomes, struct Genome *refGenome) {
    struct malnSet *malnSet;
    AllocVar(malnSet);
    malnSet->genomes = genomes;
    malnSet->refGenome = refGenome;
    return malnSet;
}

/* Construct a malnSet from a MAF file. defaultBranchLength is used to
 * assign branch lengths when inferring trees from pair-wise MAFs. */
struct malnSet *malnSet_constructFromMaf(struct Genomes *genomes, struct Genome *refGenome, char *mafFileName, double defaultBranchLength) {
    struct malnSet *malnSet = malnSet_construct(genomes, refGenome);
    struct mafFile *mafFile = mafOpen(mafFileName);
    struct mafAli *ali;
    while ((ali = mafNext(mafFile)) != NULL) {
        malnSet_addBlk(malnSet, mafAliToMAlnBlk(genomes, refGenome, ali, defaultBranchLength));
        mafAliFree(&ali);
    }
    slReverse(&malnSet->blks);
    return malnSet;
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

/* write a block to a MAF */
static void writeBlkToMaf(struct malnBlk *blk, FILE *mafFh) {
    struct mafAli *ma = malnAliToMafAli(blk);
    mafWrite(mafFh, ma);
    mafAliFree(&ma);
}

/* write a malnSet to a MAF file  */
void malnSet_writeMaf(struct malnSet *malnSet, char *mafFileName) {
    FILE *mafFh = mustOpen(mafFileName, "w");
    mafWriteStart(mafFh, NULL);
    for (struct malnBlk *blk = malnSet->blks; blk != NULL; blk = blk->next) {
        writeBlkToMaf(blk, mafFh);
    }
    mafWriteEnd(mafFh);
    carefulClose(&mafFh);
}

/* get list of blocks.  Don't modify this list. */
struct malnBlk *malnSet_getBlocks(struct malnSet *malnSet) {
    return malnSet->blks;
}

/* add block to range map for all reference genome components */
static void addRefRanges(struct malnSet *malnSet, struct malnBlk *blk) {
    for (struct malnComp *mc = blk->comps; mc != NULL; mc = mc->next) {
        if (mc->seq->genome == malnSet->refGenome) {
            genomeRangeTreeAddValList(malnSet->refBlks, mc->seq->name, mc->chromStart, mc->chromEnd, slRefNew(blk));
        }
    }
}

/* build the range tree when needed */
static void buildRangeTree(struct malnSet *malnSet) {
    malnSet->refBlks = genomeRangeTreeNew();
    for (struct malnBlk *blk = malnSet->blks; blk != NULL; blk = blk->next) {
        addRefRanges(malnSet, blk);
    }
}

/* get list of slRefs to blks who's reference range overlaps the specified
 * range.  slFreeList of results when done. */
struct slRef *malnSet_getOverlapping(struct malnSet *malnSet, struct Seq *seq, int chromStart, int chromEnd) {
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
