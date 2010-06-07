/* join two mafInvert objects */
#include "common.h"
#include "mafInvertJoin.h"
#include "mafInvert.h"
#include "hash.h"

/* flat a list of columns are used for join */
static void flatColumnList(struct MafInvertCol *srcCols) {
    for (struct MafInvertCol *srcCol = srcCols; srcCol != NULL; srcCol = srcCol->next) {
        srcCol->joinCol = TRUE;
    }
}

/* flag join columns, which is used to know where to stop when joining indel columns */
static void flagRefSeqCommon(struct MafInvertSeq *mi1Seq, struct MafInvertSeq *mi2Seq) {
    for (int pos = 0; pos < mi1Seq->seq->size; pos++) {
        if ((mi1Seq->cols[pos] != NULL) && (mi2Seq->cols[pos] != NULL)) {
            flatColumnList(mi1Seq->cols[pos]);
            flatColumnList(mi2Seq->cols[pos]);
        }
    }
}

/* Tests if a cell should be added to the current joined column. */
static bool shouldAddCell(struct Genome *ref, 
                          struct MafInvertCol *joinCol,
                          struct MafInvertCell *miCell) {
    for (int iRow = 0; iRow < joinCol->numRows; iRow++) {
        if (mafInvertCellEq(mafInvertColGetCell(joinCol, iRow), miCell)) {
            return FALSE;
        }
    }
    return TRUE;
}

/* sum number of rows in a list of columns  */
static int mafInvertColTotalSize(struct MafInvertCol *cols) {
    int sz = 0;
    for (struct MafInvertCol *col = cols; col != NULL; col = col->next) {
        sz += col->numRows;
    }
    return sz;
}

/* copy a cell to a new column, reference column map */
static void cpCell(struct MafInvert *miJoin, struct MafInvertCol *joinCol, struct MafInvertCell *srcCell) {
    int iJoinRow = joinCol->numRows++;
    struct MafInvertCell *joinCell = mafInvertColGetCell(joinCol, iJoinRow);
    *joinCell = *srcCell;
    if (mafInvertCellIsRefAligned(miJoin->ref, joinCell)) {
        mafInvertColAddToMap(miJoin, joinCol, joinCell);
    }
    mafInvertCellAssert(joinCell);
}

/* join left-adjacent indel columns */
static void joinLeftAdjIndel(struct MafInvertCol *joinCol, struct MafInvertCol *srcCol) {
    srcCol = srcCol->leftAdj;
    while ((srcCol != NULL) && (!srcCol->done) && (!srcCol->joinCol)) {
        struct MafInvertCol *joinColIns = mafInvertColClone(srcCol);
        mafInvertColInsertLeft(joinCol, joinColIns);
        srcCol->done = TRUE;
        srcCol = srcCol->leftAdj;
        joinCol = joinColIns;
    }
}

/* join right-adjacent indel columns. */
static void joinRightAdjIndel(struct MafInvertCol *joinCol, struct MafInvertCol *srcCol) {
    srcCol = srcCol->rightAdj;
    while ((srcCol != NULL) && (!srcCol->done) && (!srcCol->joinCol)) {
        struct MafInvertCol *joinColIns = mafInvertColClone(srcCol);
        mafInvertColInsertRight(joinCol, joinColIns);
        srcCol->done = TRUE;
        srcCol = srcCol->rightAdj;
        joinCol = joinColIns;
    }
}

/* add a column to the joined maf  */
static void joinColumn(struct MafInvert *miJoin, struct MafInvertCol *joinCol, struct MafInvertCol *srcCol) {
    for (int iRow = 0; iRow < srcCol->numRows; iRow++) {
        struct MafInvertCell *srcCell = mafInvertColGetCell(srcCol, iRow);
        if (shouldAddCell(miJoin->ref, joinCol, srcCell)) {
            cpCell(miJoin, joinCol, srcCell);
        }
    }
    srcCol->done = TRUE;
    joinLeftAdjIndel(joinCol, srcCol);
    joinRightAdjIndel(joinCol, srcCol);
    mafInvertColAssert(joinCol);
}

/* add a list of columns to the joined maf */
static void joinColumnList(struct MafInvert *miJoin, struct MafInvertCol *joinCol, struct MafInvertCol *srcCols) {
    for (struct MafInvertCol *srcCol = srcCols; srcCol != NULL; srcCol = srcCol->next) {
        if (!srcCol->done) {
            joinColumn(miJoin, joinCol, srcCol);
        }
    }
}

/* get the tree for a column list, which must be the same
 * (damn, I think this is gonna break things) */
static stMafTree *getColumnListMTree(struct MafInvertCol *cols) {
    if (cols == NULL) {
        return NULL;
    }
    stMafTree *mTree = cols->mTree;
    if (mTree == NULL) {
        errAbort("NULL MAF tree in column");
    }
    for (struct MafInvertCol *col = cols; col != NULL; col = col->next) {
        if (col->mTree != mTree) {
            errAbort("not all MAF trees the same in a column list");
        }
    }
    return mTree;
}

/* Join a column of the invert maf */
static struct MafInvertCol *joinSeqColumn(struct MafInvert *miJoin, struct MafInvertCol *cols1, struct MafInvertCol *cols2, struct MafInvertCol *prevJoinCol) {
    getColumnListMTree(cols1);
    getColumnListMTree(cols2);
    struct MafInvertCol *joinCol = mafInvertColNew(mafInvertColTotalSize(cols1) + mafInvertColTotalSize(cols2), NULL);
    joinColumnList(miJoin, joinCol, cols1);
    joinColumnList(miJoin, joinCol, cols2);
    if (prevJoinCol != NULL) {
        mafInvertColInsertLeft(mafInvertColFindLeftMost(joinCol), mafInvertColFindRightMost(prevJoinCol));
    }
    return joinCol;
}

/* join columns when the reference sequences containing both columns */
static void joinRefSeqCommon(struct Seq *seq, struct MafInvert *miJoin, struct MafInvertSeq *mi1Seq, struct MafInvertSeq *mi2Seq) {
    struct MafInvertCol *prevJoinCol = NULL;
    for (int pos = 0; pos < seq->size; pos++) {
        if ((mi1Seq->cols[pos] != NULL) && (mi2Seq->cols[pos] != NULL)) {
            prevJoinCol = joinSeqColumn(miJoin, mi1Seq->cols[pos], mi2Seq->cols[pos], prevJoinCol);
        }
    }
}

/* Join a sequence of the reference genome */
static void joinRefSeq(struct Seq *seq, struct MafInvert *miJoin, struct MafInvert *mi1, struct MafInvert *mi2) {
    // either of these can come back NULL, since seq may not be in both
    // alignments from the MAF
    struct MafInvertSeq *mi1Seq = hashFindVal(mi1->refSeqMap, seq->name);
    struct MafInvertSeq *mi2Seq = hashFindVal(mi2->refSeqMap, seq->name);
    assert((mi1Seq != NULL) || (mi2Seq != NULL));
    if ((mi1Seq != NULL) && (mi2Seq != NULL)) {
        flagRefSeqCommon(mi1Seq, mi2Seq);
        joinRefSeqCommon(seq, miJoin, mi1Seq, mi2Seq);
    }
}

/* Join two inverted mafs based on the reference sequence */
struct MafInvert *mafInvertJoin(struct Genomes *genomes, struct MafInvert *mi1, struct MafInvert *mi2) {
    assert(mi1->ref == mi2->ref);
    struct MafInvert *miJoin = mafInvertNew(mi1->ref);
    for (struct Seq *seq = miJoin->ref->seqs; seq != NULL; seq = seq->next) {
        joinRefSeq(seq, miJoin, mi1, mi2);
    }
    return miJoin;
}
