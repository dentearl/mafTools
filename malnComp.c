#include "malnComp.h"
#include "malnCompCursor.h"
#include "malnBlk.h"
#include "mafTree.h"
#include "genome.h"
#include "common.h"
#include "stSafeC.h"
#include "dystring.h"
#include "dnautil.h"
#include "sonLibString.h"

/* get tree location for component */
enum mafTreeLoc malnComp_getLoc(struct malnComp *comp) {
    return mafTreeNodeCompLink_getLoc(comp->ncLink);
}

/* can two components be joined? */
bool malnComp_canJoin(struct malnComp *comp1, struct malnComp *comp2) {
    // need to be overlapping or adjacent and one component must be a root to join
    return malnComp_overlapAdjacent(comp1, comp2) && ((malnComp_getLoc(comp1) == mafTreeLocRoot) || (malnComp_getLoc(comp2) == mafTreeLocRoot));
}

/* compare two components to see if they overlap or are adjacent when the second 
 * component is in the specified relative orientation (-1 == reverse complemented) */
bool malnComp_overlapAdjacentOrient(struct malnComp *comp, struct malnComp *comp2, int orient) {
    int start2 = comp2->start, end2 = comp2->end;
    if (orient < 0) {
        reverseIntRange(&start2, &end2, comp2->seq->size);
    }
    return (comp->seq == comp2->seq) && (comp->start <= end2) && (comp->end >= start2);
}


/* count aligned positions */
int malnComp_countAligned(const struct malnComp *comp) {
    int c = 0;
    for (const char *p = comp->alnStr->string; *p != '\0'; p++) {
        if (isBase(*p)) {
            c++;
        }
    }
    return c;
}

/* basic component constructor using all or part of an alignment string . */
static struct malnComp *malnComp_make(struct Seq *seq, char strand, int start, int end, char *alnStr, int strStart, int strEnd) {
    assert(end-start <= strEnd-strStart);
    struct malnComp *comp;
    AllocVar(comp);
    comp->seq = seq;
    comp->strand = strand;
    comp->start = start;
    comp->end = end;
    comp->chromStart = start;
    comp->chromEnd = end;
    if (strand == '-') {
        reverseIntRange(&comp->chromStart, &comp->chromEnd, seq->size);
    }
    comp->alnStr = dyStringNew(strEnd-strStart);
    dyStringAppendN(comp->alnStr, alnStr+strStart, strEnd-strStart);

    // sanity check, as bad sequence can really break things
    if (malnComp_countAligned(comp) != (comp->end - comp->start)) {
        malnComp_dump(comp, stderr, "bases in component (%d) doesn't match range (%d)",malnComp_countAligned(comp) , end-start);
        errAbort("invalid MAF");
    }
    return comp;
}

/* component constructor */
struct malnComp *malnComp_construct(struct Seq *seq, char strand, int start, int end, char *alnStr) {
    return malnComp_make(seq, strand, start, end, alnStr, 0, strlen(alnStr));
}

/* component constructor to clone another component */
struct malnComp *malnComp_constructClone(struct malnComp *srcComp) {
    return malnComp_construct(srcComp->seq, srcComp->strand, srcComp->start, srcComp->end, malnComp_getAln(srcComp));
}

/* destructor */
void malnComp_destruct(struct malnComp *comp) {
    if (comp != NULL) {
        if (comp->ncLink != NULL) {
            comp->ncLink->comp = NULL;
        }
        dyStringFree(&comp->alnStr);
        freeMem(comp);
    }
}

/* Clear sequence to reduce memory associated with dying blocks */
void malnComp_freeSeqMem(struct malnComp *comp) {
    dyStringFree(&comp->alnStr);
}

/* component reverse complement */
struct malnComp *malnComp_reverseComplement(struct malnComp *comp) {
    int start = comp->start, end = comp->end;
    reverseIntRange(&start, &end, comp->seq->size);
    struct malnComp *rc = malnComp_construct(comp->seq, ((comp->strand == '-') ? '+' : '-'), start, end, malnComp_getAln(comp));
    reverseComplement(rc->alnStr->string, rc->alnStr->stringSize);
    return rc;
}

/* convert a chrom range to a strand range for this component */
void malnComp_chromRangeToStrandRange(struct malnComp *comp, int chromStart, int chromEnd, int *start, int *end) {
    assert(chromStart <= chromEnd);
    assert((comp->chromStart <= chromStart) && (chromEnd <= comp->chromEnd));
    *start = chromStart;
    *end = chromEnd;
    if (comp->strand == '-') {
        reverseIntRange(start, end, comp->seq->size);
    }
}

/* convert an alignment range to a sequence range, set range to 0-0 and return
 * false if none aligned */
bool malnComp_alnRangeToSeqRange(struct malnComp *comp, int alnStart, int alnEnd, int *start, int *end) {
    struct malnCompCursor cc = malnCompCursor_make(comp);
    malnCompCursor_setAlignCol(&cc, alnStart);

    // find start
    if (!malnCompCursor_isAligned(&cc)) {
        if (!malnCompCursor_advanceToNextPos(&cc)) {
            *start = *end = 0;
            return false;
        }
    }
    *start = cc.pos;
    
    // find end
    malnCompCursor_setAlignCol(&cc,  alnEnd-1);
    *end = cc.pos;
    return true;
}

/* convert a sequence range to an alignment range, set range to 0-0 and return
 * false if none aligned */
bool malnComp_seqRangeToAlnRange(struct malnComp *comp, int start, int end, int *alnStart, int *alnEnd) {
    struct malnCompCursor cc = malnCompCursor_make(comp);
    malnCompCursor_setSeqPos(&cc, start);

    // find start
    if (!malnCompCursor_isAligned(&cc)) {
        if (!malnCompCursor_advanceToNextPos(&cc)) {
            *alnStart = *alnEnd = 0;
            return false;
        }
    }
    *alnStart = cc.alnIdx;
    
    // find end
    malnCompCursor_setSeqPos(&cc,  end);
    *alnEnd = cc.alnIdx;
    return true;
}

/* convert a sequence chrom range to an alignment range, set range to 0-0 and
 * return false if none aligned */
bool malnComp_seqChromRangeToAlnRange(struct malnComp *comp, int chromStart, int chromEnd, int *alnStart, int *alnEnd) {
    int start = chromStart, end = chromEnd;
    if (comp->strand == '-') {
        reverseIntRange(&start, &end, comp->seq->size);
    }
    bool anyAligned = malnComp_seqRangeToAlnRange(comp, start, end, alnStart, alnEnd);
    if (anyAligned && (comp->strand == '-')) {
        reverseIntRange(alnStart, alnEnd, comp->alnStr->stringSize);
    }
    return anyAligned;
}

/* pad component to the specified alignment width */
void malnComp_pad(struct malnComp *comp, int width) {
    assert(malnComp_getWidth(comp) <= width);
    dyStringAppendMultiC(comp->alnStr, '-', width-malnComp_getWidth(comp));
}

/* append a column */
void malnComp_appendCol(struct malnComp *comp, char src) {
    int nbases = isBase(src) ? 1 : 0;
    dyStringAppendC(comp->alnStr, src);
    comp->end += nbases;
    if (comp->strand == '+') {
        comp->chromEnd += nbases;
    } else {
        comp->chromStart -= nbases;
    }
}


/* append the specified number of characters from the string, adjusting component
 * coordinates */
void malnComp_append(struct malnComp *comp, char *src, int len) {
    int nbases = 0;
    for (int i = 0; i < len; i++) {
        if (isBase(src[i])) {
            nbases++;
        }
    }
    dyStringAppendN(comp->alnStr, src, len);
    comp->end += nbases;
    if (comp->strand == '+') {
        comp->chromEnd += nbases;
    } else {
        comp->chromStart -= nbases;
    }
}

/* assert that a range is being extended */
static void assertExtendRange(struct malnComp *comp, struct malnComp *srcComp, int alnStart, int alnEnd) {
    // this is very expensive
#if defined(ASSERT_SLOW) && !defined(NDEBUG)
    int start, end;
    if (malnComp_alnRangeToSeqRange(srcComp, alnStart, alnEnd, &start, &end)) {
        assert(start == comp->end);
    }
#endif
}

/* Append the specified alignment region of a component to this component.
 * The component must extend the range of this component. */
void malnComp_appendCompAln(struct malnComp *comp, struct malnComp *srcComp, int alnStart, int alnEnd) {
    assert(srcComp->seq == comp->seq);
    assert(srcComp->strand == comp->strand);
    assertExtendRange(comp, srcComp, alnStart, alnEnd);
    malnComp_append(comp, malnComp_getAln(srcComp)+alnStart, alnEnd-alnStart);
}

/* Append the current column from a cursor. */
void malnComp_appendColCursor(struct malnComp *comp, struct malnCompCursor *srcCursor) {
    assert(srcCursor->comp->seq == comp->seq);
    assert(srcCursor->comp->strand == comp->strand);
    assert(srcCursor->pos == comp->end);   
    malnComp_appendCol(comp, malnCompCursor_getCol(srcCursor));
}

/* assert that the component is set-consistent */
void malnComp_assert(struct malnComp *comp) {
#ifndef NDEBUG
    assert((comp->strand == '+') || (comp->strand == '-'));
    if (comp->strand == '+') {
        assert(comp->start == comp->chromStart);
        assert(comp->end == comp->chromEnd);
    } else {
        assert(comp->start == (comp->seq->size - comp->chromEnd));
        assert(comp->end == (comp->seq->size - comp->chromStart));
    }
    assert(comp->alnStr->stringSize == comp->blk->alnWidth);
    if (malnComp_countAligned(comp) != (comp->end - comp->start)) { // FIXME
        malnBlk_dump(comp->blk, stderr, "countAligned(comp) != (comp->end - comp->start)");
    }
    assert(malnComp_countAligned(comp) == (comp->end - comp->start));
    if (comp->ncLink == NULL) {
        malnBlk_dump(comp->blk, stderr, "NULL ncLink");
    }
    assert(comp->ncLink != NULL);
    assert(comp->ncLink->comp == comp);
    mafTreeNodeCompLink_assert(comp->ncLink);
#endif
}

/* compare two components for deterministic sorting */
int malnComp_cmp(struct malnComp *comp1, struct malnComp *comp2) {
    int diff = seqCmp(comp1->seq, comp2->seq);
    if (diff == 0) {
        diff = comp1->strand - comp2->strand;
        if (diff == 0) {
            diff = comp1->start - comp2->start;
        }
        if (diff == 0) {
            diff = comp1->end - comp2->end;
        }
    }
    return diff;
}

/* compare two components in chrom order */
int malnComp_chromCmp(struct malnComp *comp1, struct malnComp *comp2) {
    int diff = seqCmp(comp1->seq, comp2->seq);
    if (diff == 0) {
        diff = comp1->chromStart - comp2->chromStart;
    }
    if (diff == 0) {
        diff = comp1->chromEnd - comp2->chromEnd;
    }
    return diff;
}

/* construct an component from a subrange of this component. Return NULL if
 * the subrange does not contain any aligned bases. */
struct malnComp *malnComp_constructSubrange(struct malnComp *comp, int alnStart, int alnEnd) {
    struct malnCompCursor cursor = malnCompCursor_make(comp);
    malnCompCursor_setAlignCol(&cursor, alnStart);
    int compStart = cursor.pos;
    malnCompCursor_setAlignCol(&cursor, alnEnd);
    int compEnd = cursor.pos;
    struct malnComp *subComp = NULL;
    if (compStart < compEnd) {
        subComp = malnComp_make(comp->seq, comp->strand, compStart, compEnd, comp->alnStr->string, alnStart, alnEnd);
    }
    return subComp;
}

/* count the number of aligned bases in a range */
int malnComp_countAlignedRange(struct malnComp *comp, int alnStart, int alnEnd) {
    assert((0 <= alnStart) && (alnStart < alnEnd) && (alnEnd <= comp->blk->alnWidth));

    int cnt = 0;
    for (int iCol = alnStart; iCol < alnEnd; iCol++) {
        if (isBase(comp->alnStr->string[iCol])) {
            cnt++;
        }
    }
    return cnt;
}

/* are there any aligned bases in a range */
bool malnComp_anyAlignedRange(struct malnComp *comp, int alnStart, int alnEnd) {
    assert((0 <= alnStart) && (alnStart < alnEnd) && (alnEnd <= comp->blk->alnWidth));

    for (int iCol = alnStart; iCol < alnEnd; iCol++) {
        if (isBase(comp->alnStr->string[iCol])) {
            return true;
        }
    }
    return false;
}

/* find the previous aligned column, at or before a starting point, return
 * -1 if no more */
int malnComp_findPrevAligned(struct malnComp *comp, int alnStart) {
    assert((0 <= alnStart) && (alnStart < comp->blk->alnWidth));

    int iCol = alnStart;
    for (; iCol >= 0; iCol--) {
        if (isBase(comp->alnStr->string[iCol])) {
            break;
        }
    }
    return iCol;
}

/* find the next aligned column, at or after a starting point, return
 * alnWidth if no more */
int malnComp_findNextAligned(struct malnComp *comp, int alnStart) {
    assert((0 <= alnStart) && (alnStart < comp->blk->alnWidth));

    int iCol = alnStart;
    for (; iCol < comp->blk->alnWidth; iCol++) {
        if (isBase(comp->alnStr->string[iCol])) {
            break;
        }
    }
    return iCol;
}

/* print base information describing a comp, newline not included */
void malnComp_prInfo(struct malnComp *comp, FILE *fh) {
    const char *loc = (comp->ncLink != NULL) ? mafTreeLoc_str(malnComp_getLoc(comp)) : "NULL";
    fprintf(fh, "%s (%c) %d-%d %d-%d [%d] %s", comp->seq->orgSeqName, comp->strand, comp->start, comp->end, 
            ((comp->strand == '-') ? comp->chromEnd : comp->chromStart),
            ((comp->strand == '-') ? comp->chromStart : comp->chromEnd),
            malnComp_getAligned(comp), loc);
}

/* print a component for debugging purposes */
void malnComp_dumpv(struct malnComp *comp, FILE *fh, const char *label, va_list args) {
    char *fmtLabel = stSafeCDynFmtv(label, args);
    fprintf(fh, "%s ", fmtLabel);
    malnComp_prInfo(comp, fh);
    fprintf(fh, " %s\n", comp->alnStr->string);
    freeMem(fmtLabel);
}

/* print a component for debugging purposes */
void malnComp_dump(struct malnComp *comp, FILE *fh, const char *label, ...) {
    va_list args;
    va_start(args, label);
    malnComp_dumpv(comp, fh, label, args);
    va_end(args);
}
