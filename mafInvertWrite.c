/* convert a mafInvert object back to a MAF file */
#include "common.h"
#include "mafInvertWrite.h"
#include "mafInvert.h"
#include "stMafTree.h"
#include "sonLibETree.h"
#include "genome.h"
#include "jkmaf.h"
#include "dystring.h"

/* a current component in the block being built */
struct BlkComp {
    struct BlkComp *next;
    struct Seq *seq;        /* sequence currently being assembled. */
    char strand;
    int startPos;           /* -1 if not set */
    int nextPos;
    struct dyString *text;
};

/* current block being assembled */
struct BlkBuilder {
    struct Genome *refGenome;
    stMafTree *mTree;
    int width;             /* current number of columns */
    struct BlkComp *comps; /* components being built */
};

/* constructor */
static struct BlkComp *blkCompNew(struct Seq *seq, char strand, int startPos) {
    struct BlkComp *bc;
    AllocVar(bc);
    bc->seq = seq;
    bc->strand = strand;
    bc->startPos = startPos;
    bc->nextPos = startPos;
    bc->text = dyStringNew(256);
    return bc;
}

/* destructor */
static void blkCompFree(struct BlkComp *bc) {
    dyStringFree(&bc->text);
    freeMem(bc);
}

#if 0 // FIXME:
/* get length of a blkComp? */
static int blkCompLen(const struct BlkComp *bc) {
    return bc->nextPos - bc->startPos;
}
#endif

/* can a blkComp be extended to the specified position? */
static bool blkCompCanExtend(const struct BlkComp *bc, const struct Seq *seq, char strand, int pos) {
    return (bc->seq == seq) && (bc->strand == strand) && (bc->nextPos == pos);
}

/* fill a blkComp with del chars to the specified width */
static void blkCompFill(struct BlkComp *bc, int width) {
    assert(width >= bc->text->stringSize);
    dyStringAppendMultiC(bc->text, '-', width - bc->text->stringSize);
}

/* append a cell to a blkComp */
static void blkCompAppend(struct BlkComp *bc, struct MafInvertCell *cell) {
    dyStringAppendC(bc->text, cell->base);
    if (isBase(cell->base)) {
        assert(bc->nextPos == cell->pos);
        bc->nextPos++;
    }
}

/* convert a blkComp to a mafComp */
static struct mafComp *blkCompToMafComp(struct BlkComp *bc) {
    struct mafComp *mc;
    AllocVar(mc);
    mc->src = seqMkName(bc->seq);
    mc->srcSize = bc->seq->size;
    mc->strand = bc->strand;
    mc->start = bc->startPos;
    mc->size = bc->nextPos - bc->startPos;
    mc->text = cloneString(bc->text->string);
    return mc;
}

/* constructor */
static struct BlkBuilder *blkBuilderNew(struct Genome *refGenome) {
    struct BlkBuilder *bb;
    AllocVar(bb);
    bb->refGenome = refGenome;
    return bb;
}

/* destructor */
static void blkBuilderFree(struct BlkBuilder *bb) {
    struct BlkComp *bc;
    while ((bc = slPopHead(&bb->comps)) != NULL) {
        blkCompFree(bc);
    }
    freeMem(bb);
}

/* search for a BlkComp object that will extended position */
static struct BlkComp *blkBuilderGetComp(struct BlkBuilder *bb, struct Seq *seq, char strand, int pos) {
    for (struct BlkComp *bc = bb->comps; bc != NULL; bc = bc->next) {
        if (blkCompCanExtend(bc, seq, strand, pos)) {
            return bc;
        }
    }
    return NULL;
}

/* Obtain a BlkComp object that will extended position, creating and padding to
 * the current alignment length if it doesn't exist */
static struct BlkComp *blkBuilderObtainComp(struct BlkBuilder *bb, struct Seq *seq, char strand, int pos) {
    struct BlkComp *bc = blkBuilderGetComp(bb, seq, strand, pos);
    if (bc == NULL) {
        bc = blkCompNew(seq, strand, pos);
        slAddHead(&bb->comps, bc);
        blkCompFill(bc, bb->width);
    }
    return bc;
}

/* Add a cell to the block */
static void blkBuilderAddCell(struct BlkBuilder *bb, struct MafInvertCell *cell) {
    struct BlkComp *bc = blkBuilderObtainComp(bb, cell->seq, cell->strand, cell->pos);
    blkCompAppend(bc, cell);
}

/* pad alignments out to the current width only added at most one del */
static void blkBuilderPadBy1(struct BlkBuilder *bb) {
    for (struct BlkComp *bc = bb->comps; bc != NULL; bc = bc->next) {
        assert((bc->text->stringSize == bb->width) || (bc->text->stringSize == bb->width-1));
        if (bc->text->stringSize < bb->width) {
            dyStringAppendC(bc->text, '-');
        }
    }
}

/* Add a column to the block, assumes all linked blocks are added at once */
static void blkBuilderAddCol(struct BlkBuilder *bb, struct MafInvertCol *col) {
    assert(!col->done);
    for (int iRow = 0; iRow < col->numRows; iRow++) {
        blkBuilderAddCell(bb, mafInvertColGetCell(col, iRow));
    }
    bb->width++;
    blkBuilderPadBy1(bb);
    if (col->mTree != NULL) {
        assert((bb->mTree == NULL) || (bb->mTree == col->mTree));
        bb->mTree = col->mTree;
    }
    col->done = TRUE;
}

/* Fill in a BlkBuilder, start with the block containing the specified column,
 * including adjacencies not in the ref seq array. */
static struct BlkBuilder *blkBuilderFillFromRefSeq(struct MafInvertSeq *refSeq, int pos) {
    struct BlkBuilder *bb = blkBuilderNew(refSeq->seq->genome);
    for (struct MafInvertCol *col = mafInvertColFindLeftMost(refSeq->cols[pos]); col != NULL; col = col->rightAdj) {
        blkBuilderAddCol(bb, col);
    }
    return bb;
}

/* globals for compare */
static stMafTree *gMafTree;

/* compare function for two BlkComp objects by tree order */
static int blkCompCmp(const void *vbc1, const void *vbc2) {
    const struct BlkComp *bc1 = *((const struct BlkComp**)vbc1);
    const struct BlkComp *bc2 = *((const struct BlkComp**)vbc2);
    char *seq1Name = seqMkName(bc1->seq);
    char *seq2Name = seqMkName(bc2->seq);
    int order1 = stMafTree_getMafOrder(gMafTree, seq1Name);
    int order2 = stMafTree_getMafOrder(gMafTree, seq2Name);
    free(seq1Name);
    free(seq2Name);
    return order1 - order2;
}

/* sort components, putting reference sequences first */
static void blkBuilderSort(struct BlkBuilder *bb) {
    assert(bb->mTree != NULL);
    gMafTree = bb->mTree;
    slSort(&bb->comps, blkCompCmp);
}

/* convert a BlkBuilder to a MafAli */
static struct mafAli *blkBuilderToMafAli(struct BlkBuilder *bb) {
    struct mafAli *ma;
    AllocVar(ma);
    if (bb->mTree != NULL) {
        ma->tree = eTree_getNewickTreeString(stMafTree_getETree(bb->mTree));
    }
    for (struct BlkComp *bc = bb->comps; bc != NULL; bc = bc->next) {
        if (bc->startPos >= 0) {
            struct mafComp *mc = blkCompToMafComp(bc);
            if (ma->components == NULL) {
                ma->textSize = strlen(mc->text);
            }
            slAddHead(&ma->components, mc);
            assert(strlen(mc->text) == ma->textSize);
        }
    }
    slReverse(&ma->components);
    return ma;
}

/* build and output one maf block */
static void blkBuilderBuildBlk(struct MafInvertSeq *refSeq, int pos, FILE *mafFh) {
    struct BlkBuilder *bb = blkBuilderFillFromRefSeq(refSeq, pos);
    blkBuilderSort(bb);
    struct mafAli *ma = blkBuilderToMafAli(bb);
    mafWrite(mafFh, ma);
    mafAliFree(&ma);
    blkBuilderFree(bb);
}

/* Skip columns with no data or marked as done */
static int blkBuilderSkipToPendCol(struct MafInvertSeq *refSeq, int pos) {
    while ((pos < refSeq->seq->size) && ((refSeq->cols[pos] == NULL) || refSeq->cols[pos]->done)) {
        pos++;
    }
    return pos;
}

/* build and output maf blocks for one sequence in the reference genome */
static void blkBuilderBuildSeq(struct MafInvertSeq *refSeq, FILE *mafFh) {
    int pos = 0;
    while ((pos = blkBuilderSkipToPendCol(refSeq, pos)) < refSeq->seq->size) {
        blkBuilderBuildBlk(refSeq, pos, mafFh);
        pos++;
    }
}

/* convert a MafInvert to a MAF */
void mafInvertToMaf(struct MafInvert *mi, char *mafFile) {
    FILE *mafFh = mustOpen(mafFile, "w");
    mafWriteStart(mafFh, NULL);
    for (struct MafInvertSeq *refSeq = mi->refSeqs; refSeq != NULL; refSeq = refSeq->next) {
        blkBuilderBuildSeq(refSeq, mafFh);
    }
    mafWriteEnd(mafFh);
    carefulClose(&mafFh);
}
