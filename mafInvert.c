/* inverted representation of a MAF around one reference sequence */ 
#include "common.h"
#include "mafInvert.h"
#include "hash.h"
#include "jkmaf.h"
#include "obscure.h"

/* print a cell for debugging purposes */
void mafInvertCellDump(struct MafInvertCell *cell, FILE *fh) {
    fprintf(fh, "\t%s.%s %c %d %c\n", cell->seq->genome->name, cell->seq->name,
            cell->strand, cell->pos, cell->base);
}

/* constructor */
struct MafInvertCol *mafInvertColNew(int maxNumRows) {
    struct MafInvertCol *col;
    AllocVar(col);
    col->cells = needMem(maxNumRows * sizeof(struct MafInvertCell));
    return col;
}

/* clone constructor. Only clones cells, not adjacencies */
struct MafInvertCol *mafInvertColClone(struct MafInvertCol *srcCol) {
    struct MafInvertCol *col;
    AllocVar(col);
    col->numRows = srcCol->numRows;
    col->cells = needMem(col->numRows * sizeof(struct MafInvertCell));
    memcpy(col->cells, srcCol->cells, col->numRows * sizeof(struct MafInvertCell));
    return col;
}

/* print a column for debugging purposes */
void mafInvertColDump(struct MafInvertCol *col, FILE *fh, char *label) {
    fprintf(fh, "%s 0x%x\n", label, ptToInt(col));
    for (int iRow = 0; iRow < col->numRows; iRow++) {
        mafInvertCellDump(mafInvertColGetCell(col, iRow), fh);
    }
}

/* constructor */
static struct MafInvertSeq *mafInvertSeqNew(struct Seq *seq) {
    struct MafInvertSeq *miSeq;
    AllocVar(miSeq);
    miSeq->seq = seq;
    miSeq->cols = needLargeZeroedMem(seq->size * sizeof(struct MafInvertCol*));
    return miSeq;
}

/* constructor */
struct MafInvert *mafInvertNew(struct Genome *ref) {
    struct MafInvert *mi;
    AllocVar(mi);
    mi->ref = ref;
    mi->refSeqMap = hashNew(8);
    return mi;
}

/* obtain a MafInvertSeq object, creating if it doesn't exist */
static struct MafInvertSeq *mafInvertObtainSeq(struct MafInvert *mi, struct Seq *seq) {
    struct MafInvertSeq *miSeq = hashFindVal(mi->refSeqMap, seq->name);
    if (miSeq == NULL) {
        miSeq = mafInvertSeqNew(seq);
        hashAdd(mi->refSeqMap, seq->name, miSeq);
        slAddHead(&mi->refSeqs, miSeq);
    }
    return miSeq;
}

/* add a column and reference sequence cell to the ref position map */
void mafInvertColAddToMap(struct MafInvert *mi, struct MafInvertCol *col, struct  MafInvertCell *cell) {
    assert(cell->seq->genome == mi->ref);
    assert(cell->pos >= 0);
    struct MafInvertSeq *refSeq = mafInvertObtainSeq(mi, cell->seq);
    slAddHead(&(refSeq->cols[cell->pos]), col);
}

/* information about the rows being processed in a mafComp */
struct rowInfo {
    struct mafComp *comp;
    struct Seq *seq;     // genome sequence
    int nextPos;         // next position in sequence (pos strand)
    char strand;
};

/* information about block being assembled */
struct blkInfo {
    int numRows;
    struct rowInfo *rows;
    struct MafInvertCol *prevCol; // previous column reference
};

/* fill in a row info object from a mafComp */
static void rowInfoFromMafComp(struct Genomes *genomes, struct mafComp *comp, struct rowInfo *row) {
    char buf[128];
    char *srcDb = mafCompGetSrcDb(comp, buf, sizeof(buf));
    if (srcDb == NULL) {
        errAbort("Error: no org name in MAF component, source must be org.seq: %s", comp->src);
    }
    row->comp = comp;
    row->seq = genomesObtainSeq(genomes, srcDb, mafCompGetSrcName(comp), comp->srcSize);
    row->nextPos = comp->start;
    row->strand = comp->strand;
}

/* initialize blkInfo object on the stack, rows should already be allocated on the stack */
static void blkInfoInit(struct Genomes *genomes, struct mafAli *ali, int numRows, struct blkInfo *blk, struct rowInfo *rows) {
    blk->numRows = numRows;
    blk->rows = rows;
    blk->prevCol = NULL;
    struct mafComp *comp = ali->components;
    for (int iRow = 0; iRow < numRows; iRow++, comp = comp->next) {
        rowInfoFromMafComp(genomes, comp, &(rows[iRow]));
    }
}

/* initialize a MafInvertCell object from an alignment column */
static void mkMafInvertCell(struct MafInvert *mi, struct MafInvertCol *col, struct rowInfo *row, int iAli) {
    col->numRows++;
    struct MafInvertCell *cell = mafInvertColGetCell(col, col->numRows-1);
    cell->seq = row->seq;
    cell->base = row->comp->text[iAli];
    cell->pos = row->nextPos++;
    cell->strand = row->comp->strand;
    if (mafInvertCellIsRefAligned(mi->ref, cell)) {
        mafInvertColAddToMap(mi, col, cell);
    }
    assert(isBase(cell->base));
}

/* initialize cells in a MafInvertCol */
static void mkMafInvertCells(struct MafInvert *mi, struct MafInvertCol *col, struct blkInfo *blk, int iAli) {
    for (int iRow = 0; iRow < blk->numRows; iRow++) {
        struct rowInfo *row = &(blk->rows[iRow]);
        if (isBase(row->comp->text[iAli])) {
            mkMafInvertCell(mi, col, row, iAli);
        }
    }
}

/* construct a MafInvertCol object from an alignment column */
static void mkMafInvertCol(struct MafInvert *mi, struct blkInfo *blk, int iAli, ETree *tree) {
    struct MafInvertCol *col = mafInvertColNew(blk->numRows);
    if (blk->prevCol != NULL) {
        mafInvertColInsertLeft(col, blk->prevCol);
    }
    mkMafInvertCells(mi, col, blk, iAli);
    blk->prevCol = col;
}

/* get the tree for a mafAli, either building it from the data
 * or parsing it. */
static ETree *getMafAliTree(struct MafInvert *mi, struct mafAli *ali) {
    if (ali->tree == NULL) {
        errAbort("mafAli does not have a tree"); // FIXME: tmp
        return NULL;
    } else {
        return eTree_parseNewickString(ali->tree);
    }
}

/* Add columns from a mafAli to the inverted maf */
static void mafInvertAddMafAli(struct MafInvert *mi, struct Genomes *genomes, struct mafAli *ali) {
    ETree *tree = getMafAliTree(mi, ali);
    int numRows = slCount(ali->components);
    struct blkInfo blk;
    struct rowInfo rows[numRows];
    blkInfoInit(genomes, ali, numRows, &blk, rows);
    for (int iAli = 0; iAli < ali->textSize; iAli++) {
        mkMafInvertCol(mi, &blk, iAli, tree);
    }
}

/* Build an inverted MAF from a MAF file */
struct MafInvert *mafInvertBuild(struct Genomes *genomes, struct Genome *ref, char *fileName) {
    struct mafFile *mafFile = mafOpen(fileName);
    struct MafInvert *mi = mafInvertNew(ref);
    struct mafAli *ali;
    while ((ali = mafNext(mafFile)) != NULL) {
        mafInvertAddMafAli(mi, genomes, ali);
        mafAliFree(&ali);
    }
    return mi;
}
