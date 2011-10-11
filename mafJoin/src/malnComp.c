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
#include "sonLibTree.h"

/* construct a new seg. base maybe NULL */
static struct malnCompSeg *malnCompSeg_construct(int start, int alnStart, char *bases, int len) {
    assert((start >= 0) && (alnStart >= 0) && (len >= 0) && (len <= strlen(bases)));
    struct malnCompSeg *seg;
    AllocVar(seg);
    seg->start = start;
    seg->alnStart = alnStart;
    if (bases != NULL) {
        seg->bases = cloneStringZ(bases, len);
        seg->size = len;
    }
    return seg;
}

/* clone a segment */
struct malnCompSeg *malnCompSeg_constructClone(struct malnCompSeg *srcSeg) {
    struct malnCompSeg *seg;
    AllocVar(seg);
    seg->start = srcSeg->start;
    seg->alnStart = srcSeg->alnStart;
    if (srcSeg->bases != NULL) {
        seg->bases = cloneString(srcSeg->bases);
        seg->size = srcSeg->size;
    }
    return seg;
}

/* extend a segment */
void malnCompSeg_extend(struct malnCompSeg *seg, char *newBases, int newSize) {
    seg->bases = needMoreMem(seg->bases, seg->size+1, seg->size+newSize+1);
    strcpy(seg->bases+seg->size, newBases);
    seg->size += newSize;
}

/* clone, creating a reverse complement of a segment */
static struct malnCompSeg *malnCompSeg_reverseComplement(struct malnCompSeg *srcSeg, struct Seq *seq, int alnWidth) {
    int start = seq->size - (srcSeg->start + srcSeg->size);
    int alnStart = alnWidth - (srcSeg->alnStart + srcSeg->size);
    struct malnCompSeg *seg = malnCompSeg_construct(start, alnStart, srcSeg->bases, srcSeg->size);
    reverseComplement(seg->bases, seg->size);
    return seg;
}

/* destruct a segment */
static void malnCompSeg_destruct(struct malnCompSeg *seg) {
    if (seg != NULL) {
        freeMem(seg->bases);
        freeMem(seg);
    }
}


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

/* find next contiguous range of bases */
static int nextContigBases(char *alnStr, int *alnNext, int alnEnd) {
    int i = *alnNext;
    // find next base
    while ((i < alnEnd) && !isBase(alnStr[i])) {
        i++;
    }

    int rngStart = i;

    // find end of range
    while ((i < alnEnd) && isBase(alnStr[i])) {
        i++;
    }
    *alnNext = i;
    return rngStart;
}

/* build segments from sequence */
static struct malnCompSeg *alnStrToSegs(int start, int end, char *alnStr, int alnStart, int alnEnd) {
    struct malnCompSeg *segs = NULL;
    int nextPos = start;
    int alnNext = alnStart, alnRngStart;
    while ((alnRngStart = nextContigBases(alnStr, &alnNext, alnEnd)) < alnEnd) {
        slAddHead(&segs, malnCompSeg_construct(nextPos, alnRngStart, alnStr+alnRngStart, (alnNext-alnRngStart)));
        nextPos += alnNext-alnRngStart;
    }
    if (nextPos != end) {
        errAbort("sequence range doesn't match bases in alignment");
    }
    slReverse(&segs);
    return segs;
}

/* constructor */
struct malnComp *malnComp_construct(struct Seq *seq, char strand, int start, int end, int alnWidth) {
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
    comp->alnWidth = alnWidth;
    return comp;
}

/* constructor from alignment string  */
struct malnComp *malnComp_constructFromStr(struct Seq *seq, char strand, int start, int end, char *alnStr) {
    struct malnComp *comp = malnComp_construct(seq, strand, start, end, strlen(alnStr));
    comp->segs = alnStrToSegs(start, end, alnStr, 0, comp->alnWidth);
    return comp;
}

/* component constructor to clone another component */
struct malnComp *malnComp_constructClone(struct malnComp *srcComp) {
    struct malnComp *comp = malnComp_construct(srcComp->seq, srcComp->strand, srcComp->start, srcComp->end, srcComp->alnWidth);
    for (struct malnCompSeg *srcSeg = srcComp->segs; srcSeg != NULL; srcSeg = srcSeg->next) {
        slAddHead(&comp->segs, malnCompSeg_constructClone(srcSeg));
    }
    slReverse(&comp->segs);
    return comp;
}

/* Clear sequence to reduce memory associated with dying blocks */
void malnComp_freeSeqMem(struct malnComp *comp) {
    struct malnCompSeg *seg;
    while ((seg = slPopHead(&comp->segs)) != NULL) {
        malnCompSeg_destruct(seg);
    }
}

/* destructor */
void malnComp_destruct(struct malnComp *comp) {
    if (comp != NULL) {
        if (comp->ncLink != NULL) {
            comp->ncLink->comp = NULL;
        }
        malnComp_freeSeqMem(comp);
        freeMem(comp);
    }
}

/* component reverse complement */
struct malnComp *malnComp_reverseComplement(struct malnComp *srcComp) {
    int start = srcComp->start, end = srcComp->end;
    reverseIntRange(&start, &end, srcComp->seq->size);
    struct malnComp *rcComp = malnComp_construct(srcComp->seq, ((srcComp->strand == '-') ? '+' : '-'), start, end, srcComp->alnWidth);
    for (struct malnCompSeg *srcSeg = srcComp->segs; srcSeg != NULL; srcSeg = srcSeg->next) {
        slAddHead(&rcComp->segs, malnCompSeg_reverseComplement(srcSeg, srcComp->seq, rcComp->alnWidth));
    }
    // list order was reversed by above loop
    return rcComp;
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
    if (!malnCompCursor_advanceToAligned(&cc)) {
        *start = *end = 0;
        return false;
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
    if (!malnCompCursor_advanceToAligned(&cc)) {
        *alnStart = *alnEnd = 0;
        return false;
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
        reverseIntRange(alnStart, alnEnd, comp->alnWidth);
    }
    return anyAligned;
}

/* append bases from a string */
static void appendBases(struct malnComp *comp, char *alnStr, int strStart, int numBases) {
    slAddTail(&comp->segs, malnCompSeg_construct(comp->end, comp->alnWidth, alnStr+strStart, numBases));
    comp->end += numBases;
    if (comp->strand == '+') {
        comp->chromEnd += numBases;
    } else {
        comp->chromStart -= numBases;
    }
    comp->alnWidth += numBases;
}

/* assert that a range is being extended */
static void assertExtendRange(struct malnComp *comp, struct malnComp *srcComp, int alnStart, int alnEnd) {
#ifdef ASSERT_SLOW
    assert(alnStart <= alnEnd);
    int start, end;
    if ((alnStart < alnEnd) && malnComp_alnRangeToSeqRange(srcComp, alnStart, alnEnd, &start, &end)) {
        assert(start == comp->end);
    }
#endif
}

/* extend a component with unaligned columns, advancing the cursor. */
static void appendCompGap(struct malnComp *comp, struct malnCompCursor *srcCursor, int alnEnd) {
    int extendAlnEnd = (srcCursor->seg == NULL) ? alnEnd : min(srcCursor->seg->alnStart, alnEnd);
    comp->alnWidth += extendAlnEnd - srcCursor->alnIdx;
    malnCompCursor_setAlignCol(srcCursor, extendAlnEnd);
}

/* extend a component with aligned columns, advancing the cursor. */
static void appendCompSeg(struct malnComp *comp, struct malnCompCursor *srcCursor, int alnEnd) {
    int extendAlnEnd = min(malnCompSeg_getAlnEnd(srcCursor->seg), alnEnd);
    appendBases(comp, srcCursor->seg->bases, (srcCursor->alnIdx - srcCursor->seg->alnStart), (extendAlnEnd - srcCursor->alnIdx));
    malnCompCursor_setAlignCol(srcCursor, extendAlnEnd);
}

/* Append from the current position of the specified cursor up to the
 * given alignment column, advancing the cursor   */
void malnComp_appendFromCursor(struct malnComp *comp, struct malnCompCursor *srcCursor, int alnEnd) {
    // FIXME: better function name
    assert(srcCursor->comp->seq == comp->seq);
    assert(srcCursor->comp->strand == comp->strand);
    assert(srcCursor->alnIdx <= alnEnd);
    assertExtendRange(comp, srcCursor->comp, srcCursor->alnIdx, alnEnd);
   
#ifndef NDEBUG
    int newAlnWidth = comp->alnWidth + (alnEnd - srcCursor->alnIdx);
#endif
    while (srcCursor->alnIdx < alnEnd) {
        if (!srcCursor->isAligned) {
            appendCompGap(comp, srcCursor, alnEnd);
        }
        if ((srcCursor->alnIdx < alnEnd) && srcCursor->isAligned) {
            appendCompSeg(comp, srcCursor, alnEnd);
        }
    }
    assert(comp->alnWidth == newAlnWidth);
    if (comp->blk != NULL) {
        comp->blk->alnWidth = max(comp->blk->alnWidth, comp->alnWidth);
    }
}

#ifndef NDEBUG
/* count aligned positions */
static int malnComp_countAligned(const struct malnComp *comp) {
    int c = 0;
    for (struct malnCompSeg *seg = comp->segs; seg != NULL; seg = seg->next) {
        c += seg->size;
    }
    return c;
}
#endif

/* assert that the component is set-consistent. ncLinkCheck can be false
 * for sanity check before links are constructed. */
void malnComp_assert(struct malnComp *comp, bool ncLinkCheck) {
#ifndef NDEBUG
    assert((comp->strand == '+') || (comp->strand == '-'));
    if (comp->strand == '+') {
        assert(comp->start == comp->chromStart);
        assert(comp->end == comp->chromEnd);
    } else {
        assert(comp->start == (comp->seq->size - comp->chromEnd));
        assert(comp->end == (comp->seq->size - comp->chromStart));
    }
    assert(comp->start < comp->end);
    assert(comp->end <= comp->seq->size);
    assert(comp->alnWidth >= (comp->end - comp->start));
    assert(comp->alnWidth == comp->blk->alnWidth);
    if (malnComp_countAligned(comp) != (comp->end - comp->start)) { // FIXME
        malnBlk_dump(comp->blk, stderr, malnBlk_dumpDefault, "countAligned(comp) != (comp->end - comp->start)");
    }
    assert(malnComp_countAligned(comp) == (comp->end - comp->start));
    if (ncLinkCheck) {
        if (comp->ncLink == NULL) {
            malnBlk_dump(comp->blk, stderr, malnBlk_dumpDefault, "NULL ncLink");
        }
        assert(comp->ncLink != NULL);
        assert(comp->ncLink->comp == comp);
    }

    // verify contiguous, ordered blks
    assert(comp->segs->start == comp->start);
    struct malnCompSeg *prevSeg = comp->segs;
    for (struct malnCompSeg *seg = prevSeg->next; seg != NULL; seg = seg->next) {
        assert(seg->start == malnCompSeg_getEnd(prevSeg));
        prevSeg = seg;
    }
    assert(malnCompSeg_getEnd(prevSeg) == comp->end);

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
struct malnComp *malnComp_constructSubrange(struct malnComp *srcComp, int alnStart, int alnEnd) {
    struct malnCompCursor srcCursor = malnCompCursor_make(srcComp);
    malnCompCursor_setAlignCol(&srcCursor, alnEnd);  // do end first
    int srcCompEnd = srcCursor.pos;
    malnCompCursor_setAlignCol(&srcCursor, alnStart); // leave at the start
    int srcCompStart = srcCursor.pos;
    struct malnComp *subComp = NULL;
    if (srcCompStart < srcCompEnd) {
        subComp = malnComp_construct(srcComp->seq, srcComp->strand, srcCompStart, srcCompStart, 0);
        malnComp_appendFromCursor(subComp, &srcCursor, alnEnd);
    }
    return subComp;
}

/* are there any aligned bases in a range */
bool malnComp_anyAlignedRange(struct malnComp *comp, int alnStart, int alnEnd) {
    assert((0 <= alnStart) && (alnStart < alnEnd) && (alnEnd <= comp->blk->alnWidth));

    struct malnCompCursor cursor = malnCompCursor_make(comp);
    malnCompCursor_setAlignCol(&cursor, alnStart);
    return malnCompCursor_advanceToAligned(&cursor) && (cursor.alnIdx < alnEnd);
}

/* print base information describing a comp, newline not included */
void malnComp_prInfo(struct malnComp *comp, FILE *fh) {
    const char *loc = (comp->ncLink != NULL) ? mafTreeLoc_str(malnComp_getLoc(comp)) : "NULL";
    fprintf(fh, "%s (%c) %d-%d %d-%d [%d] %s", comp->seq->orgSeqName, comp->strand, comp->start, comp->end, 
            ((comp->strand == '-') ? comp->chromEnd : comp->chromStart),
            ((comp->strand == '-') ? comp->chromStart : comp->chromEnd),
            malnComp_numAlignedBases(comp), loc);
}

/* print a component for debugging purposes */
void malnComp_dumpv(struct malnComp *comp, FILE *fh, const char *label, va_list args) {
    char *fmtLabel = stSafeCDynFmtv(label, args);
    fprintf(fh, "%s ", fmtLabel);
    freeMem(fmtLabel);
    malnComp_prInfo(comp, fh);
    // FIXME: do block at a time
    fputc(' ', fh);
    struct malnCompCursor cursor = malnCompCursor_make(comp);
    while (malnCompCursor_incr(&cursor)) {
        fputc(malnCompCursor_getCol(&cursor), fh);
    }
    fputc('\n', fh);
}

/* print a component for debugging purposes */
void malnComp_dump(struct malnComp *comp, FILE *fh, const char *label, ...) {
    va_list args;
    va_start(args, label);
    malnComp_dumpv(comp, fh, label, args);
    va_end(args);
}
