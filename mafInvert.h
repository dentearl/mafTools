/* inverted representation of a MAF around one reference sequence */ 
#ifndef mafInvert_h
#define mafInvert_h
#include "mafJoinTypes.h"
#include <ctype.h>
#include <stdbool.h>
#include <stddef.h>
#include <assert.h>
#include "genome.h"
struct mafAli;


/* one cell in one column of the MAF */
struct MafInvertCell {
    struct Seq *seq;
    int pos;
    char strand;
    char base;    // doesn't have deletion chars
};

/* one column of the MAF */
struct MafInvertCol {
    struct MafInvertCol *next;     // next in list (used in building, not adjacency)
    struct MafInvertCol *leftAdj;  // adjacent cells
    struct MafInvertCol *rightAdj;
    stMafTree *mTree;              // tree for column, shared so never freed
    int numRows;
    bool joinCol;                  // column contains a join cell (in both seqs)
    bool done;                     // have to indicate column has been consumed in some manner (hacky)
    struct MafInvertCell *cells;   // array of cells, in tree order
};

/* inverted MAF for one sequence in the reference genome */
struct MafInvertSeq {
    struct MafInvertSeq *next;
    struct Seq *seq;
    struct MafInvertCol **cols;  // array of sequence position to list of columns
                                 // containing the position.
};

/* inverted MAF, indexed by rows of the reference sequence  */
struct MafInvert {
    struct Genome *ref;
    struct hash *refSeqMap;        // map and list of sequences
    stMafTreeTbl *mafTreeTbl;      // shared MAF trees, avoid dups
    struct MafInvertSeq *refSeqs;
};

/* does a character represent a base */
static inline bool isBase(char base) {
    // n.b. isalpha doesn't return 0/1, might be out of bool range
    return isalpha(base) != 0;
}

/* are two bases the same, ignoring case */
static inline bool baseEq(char base1, char base2) {
    return toupper(base1) == toupper(base2);
}

/* assert that a cell looks sane */
static inline void mafInvertCellAssert(struct MafInvertCell *cell) {
    assert(isBase(cell->base));
    assert(cell->pos >= 0);
}

/* get a pointer to a cell in a column */
static inline struct MafInvertCell *mafInvertColGetCell(struct MafInvertCol *col, int iRow) {
    assert((0 <= iRow) && (iRow < col->numRows));
    return &(col->cells[iRow]);
}

/* are two cells the same? */
static inline bool mafInvertCellEq(struct MafInvertCell *cell1, struct MafInvertCell *cell2) {
    return (cell1->seq == cell2->seq) && (cell1->pos == cell2->pos) && (cell1->strand == cell2->strand) && baseEq(cell1->base, cell2->base);
} 

/* is a cell a reference genome cell that is aligned */
static inline bool mafInvertCellIsRefAligned(struct Genome *ref, struct MafInvertCell *cell) {
    return (cell->pos >= 0) && (cell->seq->genome == ref);
}

/* constructor */
struct MafInvertCol *mafInvertColNew(int numRows, stMafTree *mTree);

/* clone constructor. Only clones cells, not adjacencies */
struct MafInvertCol *mafInvertColClone(struct MafInvertCol *srcCol);

/* determine if a cell contains a sequence and position */
static inline bool mafInvertColContains(struct MafInvertCol *col, struct Seq *seq, char strand, int pos) {
    for (int iRow = 0; iRow < col->numRows; iRow++) {
        struct MafInvertCell *cell = mafInvertColGetCell(col, iRow);
        if ((cell->seq == seq) && (cell->strand == strand) && (cell->pos == pos)) {
            return true;
        }
    }
    return false;
}

/* find the left most adjacent column, given a starting column, returns col
 * if there are no adjacentiess. */
static inline struct MafInvertCol *mafInvertColFindLeftMost(struct MafInvertCol *col) {
    while (col->leftAdj != NULL) {
        col = col->leftAdj;
    }
    return col;
}

/* find the right most adjacent column, given a starting column, returns col
 * if there are no adjacenties. */
static inline struct MafInvertCol *mafInvertColFindRightMost(struct MafInvertCol *col) {
    while (col->rightAdj != NULL) {
        col = col->rightAdj;
    }
    return col;
}

/* find the first left column without the done flag set, returns col if there
 * are no undone adjacentiess. */
static inline struct MafInvertCol *mafInvertColFindLeftUndone(struct MafInvertCol *col) {
    while ((col->leftAdj != NULL) && (!col->leftAdj->done)) {
        col = col->leftAdj;
    }
    return col;
}

/* add insert one or more linked column as left adjacency */
static inline void mafInvertColInsertLeft(struct MafInvertCol *col, struct MafInvertCol *newCol) {
    struct MafInvertCol *newLeft = mafInvertColFindLeftMost(newCol);
    struct MafInvertCol *newRight = mafInvertColFindRightMost(newCol);
    if (col->leftAdj != NULL) {
        assert(col->leftAdj->rightAdj == col);
        col->leftAdj->rightAdj = newLeft;
        newLeft->leftAdj = col->leftAdj;
    }
    col->leftAdj = newRight;
    newRight->rightAdj = col;
}

/* add insert one or more linked column as right adjacency */
static inline void mafInvertColInsertRight(struct MafInvertCol *col, struct MafInvertCol *newCol) {
    struct MafInvertCol *newLeft = mafInvertColFindLeftMost(newCol);
    struct MafInvertCol *newRight = mafInvertColFindRightMost(newCol);
    if (col->rightAdj != NULL) {
        assert(col->rightAdj->leftAdj == col);
        col->rightAdj->leftAdj = newRight;
        newRight->rightAdj = col->rightAdj;
    }
    col->rightAdj = newLeft;
    newLeft->leftAdj = col;
}

/* does a column have the reference genome? */
static inline bool mafInvertColHasRef(struct MafInvertCol *col, struct MafInvert *mi) {
    for (int iRow = 0; iRow < col->numRows; iRow++) {
        if (col->cells[iRow].seq->genome == mi->ref) {
            return true;
        }
    }
    return false;
}

/* assert that a column looks sane */
static inline void mafInvertColAssert(struct MafInvertCol *col) {
    for (int iRow = 0; iRow < col->numRows; iRow++) {
        mafInvertCellAssert(mafInvertColGetCell(col, iRow));
    }
    assert((col->leftAdj == NULL) || (col->leftAdj->rightAdj == col));
    assert((col->rightAdj == NULL) || (col->rightAdj->leftAdj == col));
}

/* add a column and reference sequence cell to the ref position map */
void mafInvertColAddToMap(struct MafInvert *mi, struct MafInvertCol *col, struct  MafInvertCell *cell);

/* constructor */
struct MafInvert *mafInvertNew(struct Genome *ref);

/* Build an inverted MAF from a MAF file */
struct MafInvert *mafInvertBuild(struct Genomes *genomes, struct Genome *ref, char *fileName);
#endif
