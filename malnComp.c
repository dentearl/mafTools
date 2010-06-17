#include "malnComp.h"
#include "malnCompCursor.h"
#include "common.h"
#include "dnautil.h"

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
    struct malnComp *rc = malnComp_construct(comp->seq, ((comp->strand == '-') ? '+' : '-'), comp->start, comp->end, malnComp_getAln(comp));
    reverseIntRange(&rc->start, &rc->end, rc->seq->size);
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
    malnCompCursor_setSeqPos(&cc,  end-1);
    *alnEnd = cc.alnIdx+1;
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
