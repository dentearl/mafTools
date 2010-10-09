#include "malnComp.h"
#include "malnCompCursor.h"
#include "genome.h"
#include "common.h"
#include "dystring.h"
#include "dnautil.h"

/* get tree location for component */
enum mafTreeLoc malnComp_getLoc(struct malnComp *comp) {
    return mafTreeNodeCompLink_getLoc(comp->ncLink);
}

/* can a join be made at this component? */
bool malnComp_joinable(struct malnComp *comp) {
    return (malnComp_getLoc(comp) & (mafTreeLocRoot|mafTreeLocLeaf)) != 0;
}

/* can two components be joined? */
bool malnComp_canJoin(struct malnComp *comp1, struct malnComp *comp2) {
    return malnComp_overlap(comp1, comp2)
        && (((malnComp_getLoc(comp1) == mafTreeLocRoot) && (malnComp_getLoc(comp2) == mafTreeLocRoot))
            || ((malnComp_getLoc(comp1) == mafTreeLocRoot) && (malnComp_getLoc(comp2) == mafTreeLocLeaf))
            || ((malnComp_getLoc(comp1) == mafTreeLocLeaf) && (malnComp_getLoc(comp2) == mafTreeLocRoot)));
}

/* count aligned positions */
static int countAligned(const struct malnComp *comp) {
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
    int n = countAligned(comp);
    if (n != (comp->end - comp->start)) {
        errAbort("component %s %d-%d(%c): aligned positions (%d) != range (%d)", seq->orgSeqName, start, end, strand, n, end-start);
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
        dyStringFree(&comp->alnStr);
        freeMem(comp);
    }
}

/* component reverse complement */
struct malnComp *malnComp_reverseComplement(struct malnComp *comp) {
    int start = comp->start, end = comp->end;
    reverseIntRange(&start, &end, comp->seq->size);
    struct malnComp *rc = malnComp_construct(comp->seq, ((comp->strand == '-') ? '+' : '-'), start, end, malnComp_getAln(comp));
    reverseComplement(rc->alnStr->string, rc->alnStr->stringSize);
    return rc;
}

/* convert an alignment range to a sequence range, set range to 0-0 and return
 * false if none aligned */
bool malnComp_alnRangeToSeqRange(struct malnComp *comp, int alnStart, int alnEnd, int *start, int *end) {
    struct malnCompCursor cc = malnCompCursor_make(comp);
    malnCompCursor_setAlignCol(&cc, alnStart);

    // find start
    if (!malnCompCursor_isAligned(&cc)) {
        if (!malnCompCursor_nextPos(&cc)) {
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
        if (!malnCompCursor_nextPos(&cc)) {
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

/* pad component to the specified alignment width */
void malnComp_pad(struct malnComp *comp, int width) {
    assert(malnComp_getWidth(comp) <= width);
    dyStringAppendMultiC(comp->alnStr, '-', width-malnComp_getWidth(comp));
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

/* Append the specified alignment region of a component to this component.
 * The component must extend the range of this component. */
void malnComp_appendCompAln(struct malnComp *comp, struct malnComp *srcComp, int alnStart, int alnEnd) {
    malnComp_append(comp, malnComp_getAln(srcComp)+alnStart, alnEnd-alnStart);
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
    assert(countAligned(comp) == (comp->end - comp->start));
    if (comp->ncLink != NULL) {
        assert(comp->ncLink->comp == comp);
    }
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

/* print a component for debugging purposes */
void malnComp_dump(struct malnComp *comp, const char *label, FILE *fh) {
    const char *loc = (comp->ncLink != NULL) ? mafTreeLoc_str(malnComp_getLoc(comp)) : "NULL";
    fprintf(fh, "%s %s (%c) %d-%d %d-%d %s %s\n", label, comp->seq->orgSeqName, comp->strand, comp->start, comp->end, comp->chromStart, comp->chromEnd, loc, comp->alnStr->string);
}
