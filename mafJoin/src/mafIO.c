#include "mafIO.h"
#include "genome.h"
#include "malnSet.h"
#include "malnComp.h"
#include "malnBlk.h"
#include "mafTree.h"
#include "malnBlkSet.h"
#include "sonLibList.h"
#include "sonLibString.h"
#include "stSafeC.h"
#include "common.h"
#include "jkmaf.h"

/* convert a mafComp to an malnComp */
static struct malnComp *mafCompToMAlnComp(struct Genomes *genomes, struct mafComp *mafComp) {
    char buf[128];
    char *srcDb = mafCompGetSrcDb(mafComp, buf, sizeof(buf));
    if (srcDb == NULL) {
        errAbort("Error: no organism name in MAF component, sample "
                 "must be formatted as \"organism.seq\", e.g. \"hg18.chr2\". "
                 "Offending sequence name: %s", mafComp->src);
    }
    struct malnComp *comp = malnComp_construct(genomesObtainSeq(genomes, srcDb, mafCompGetSrcName(mafComp), mafComp->srcSize),
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
    return -(malnComp_getAligned(comp1) - malnComp_getAligned(comp2));
}

/* set the ordering for block were tree must be constructed */
static void orderTreeless(struct malnBlk *blk, struct Genome *treelessRootGenome) {
    sortTreelessRootGenome = treelessRootGenome;
    slSort(&(blk->comps), orderTreelessCmp);
    sortTreelessRootGenome = NULL;
}

/* convert a mafAli to an malnBlk */
struct malnBlk *mafIO_malnBlkRead(struct Genomes *genomes, struct mafAli *mafAli, double defaultBranchLength, struct Genome *treelessRootGenome, bool requireTree) {
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
    } else if (requireTree) {
        errAbort("No tree fuond in mafAli block and no treelessRoot genome option specified.");
    } else {
        slReverse(&blk->comps);
    }
    malnBlk_finish(blk);
    return blk;
}

/* Initialize the alignment starts in an array of splits. Adjust bounds of
 * partitions so that the have least some bases of the block's root aligned,
 * as require by split code */
static void splitBlockInitStarts(struct malnBlk *blk, int numParts, int partSize, int alnStarts[]) {
    struct malnComp *rootComp = malnBlk_getRootComp(blk);
    int alnNext = 0;
    int iPart = 0;
    // fill in, combining partitions if needed
    while ((iPart <= numParts) && (alnNext < blk->alnWidth)) {
        int alnEnd = min(alnNext+partSize, blk->alnWidth);
        if (malnComp_anyAlignedRange(rootComp, alnNext, alnEnd)) {
            alnStarts[iPart++] = alnNext;
        }
        alnNext = alnEnd;
    }
    int iLast = iPart-1;
    // fill in the rest with end
    for (; iPart <= numParts; iPart++) {
        alnStarts[iPart] = blk->alnWidth;
    }

    // if the last entry does have any aligned root bases, move the bound back
    if (!malnComp_anyAlignedRange(rootComp, alnStarts[iLast], alnStarts[iLast+1])) {
        alnStarts[iLast] = alnStarts[iLast+1];
    }
}

/* split a block in to roughly equal sizes pieces under the limit and add
 * them.  This is very annoying as each piece must include at least some
 * bases of the root sequence. FIXME: this is all a hack around having bogus trees */
static void splitBlock(struct malnSet *malnSet, int maxInputBlkWidth, struct malnBlk *blk) {
    int numParts = (blk->alnWidth + maxInputBlkWidth - 1) / maxInputBlkWidth;
    int partSize = (blk->alnWidth + numParts - 1) / numParts;
    int alnStarts[numParts+1];
    splitBlockInitStarts(blk, numParts, partSize, alnStarts);

    for (int iPart = 0; iPart < numParts; iPart++) {
        if (alnStarts[iPart] < alnStarts[iPart+1]) {
            malnSet_addBlk(malnSet, malnBlk_constructSubrange(blk, alnStarts[iPart], alnStarts[iPart+1]));
        }
    }
}

/* add a mafAli to the set */
static void addMafAli(struct malnSet *malnSet, struct mafAli *mafAli, int maxInputBlkWidth, double defaultBranchLength, struct Genome *treelessRootGenome) {
    struct malnBlk *blk = mafIO_malnBlkRead(malnSet_getGenomes(malnSet), mafAli, defaultBranchLength, treelessRootGenome, true);
    if (blk->alnWidth < maxInputBlkWidth) {
        malnSet_addBlk(malnSet, blk);
    } else {
        splitBlock(malnSet, maxInputBlkWidth, blk);
        malnBlk_destruct(blk);
    }
}

/* Construct a malnSet from a MAF file. defaultBranchLength is used to
 * assign branch lengths when inferring trees from the MAF. */
struct malnSet *mafIO_malnSetRead(struct Genomes *genomes, char *mafFileName, int maxInputBlkWidth, double defaultBranchLength, struct Genome *treelessRootGenome) {
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
    struct malnBlkSetIterator *iter = malnSet_getBlocks(malnSet);
    struct malnBlk *blk;
    while ((blk = malnBlkSetIterator_getNext(iter)) != NULL) {
        stList_append(sorted, blk);
    }
    malnBlkSetIterator_destruct(iter);
    stList_sort(sorted, blkCmpRootComp);
    return sorted;
}

/* write a block to a MAF */
void mafIO_malnBlkWrite(struct malnBlk *blk, FILE *mafFh) {
    malnBlk_validate(blk);
    struct mafAli *ma = malnAliToMafAli(blk);
    mafWrite(mafFh, ma);
    mafAliFree(&ma);
}

/* write a malnSet to a MAF file  */
void mafIO_malnSetWrite(struct malnSet *malnSet, char *mafFileName) {
    malnSet_assert(malnSet);
    stList *sorted = buildRootSorted(malnSet);

    FILE *mafFh = mustOpen(mafFileName, "w");
    mafWriteStart(mafFh, NULL);
    stListIterator *iter = stList_getIterator(sorted);
    struct malnBlk *blk;
    while ((blk = stList_getNext(iter)) != NULL) {
        mafIO_malnBlkWrite(blk, mafFh);
    }
    stList_destructIterator(iter);
    mafWriteEnd(mafFh);
    carefulClose(&mafFh);
    stList_destruct(sorted);
}

