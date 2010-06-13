#include "stMalnComp.h"
#include "stMalnCompCursor.h"
#include "common.h"
#include "dnautil.h"

/* basic component constructor using all or part of an alignment string . */
static struct stMalnComp *stMalnComp_make(struct Seq *seq, char strand, int start, int end, char *alnStr, int strStart, int strEnd) {
    assert(end-start <= strEnd-strStart);
    struct stMalnComp *comp;
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
struct stMalnComp *stMalnComp_construct(struct Seq *seq, char strand, int start, int end, char *alnStr) {
    return stMalnComp_make(seq, strand, start, end, alnStr, 0, strlen(alnStr));
}

/* component constructor to clone another component */
struct stMalnComp *stMalnComp_constructClone(struct stMalnComp *srcComp) {
    return stMalnComp_construct(srcComp->seq, srcComp->strand, srcComp->start, srcComp->end, stMalnComp_getAln(srcComp));
}

/* destructor */
void stMalnComp_destruct(struct stMalnComp *comp) {
    if (comp != NULL) {
        dyStringFree(&comp->alnStr);
        freeMem(comp);
    }
}

/* component reverse complement */
struct stMalnComp *stMalnComp_reverseComplement(struct stMalnComp *comp) {
    struct stMalnComp *rc = stMalnComp_construct(comp->seq, ((comp->strand == '-') ? '+' : '-'), comp->start, comp->end, stMalnComp_getAln(comp));
    reverseIntRange(&rc->start, &rc->end, rc->seq->size);
    reverseComplement(rc->alnStr->string, rc->alnStr->stringSize);
    return rc;
}

/* convert an alignment range to a sequence range, set range to 0-0 and return
 * false if none aligned */
bool stMalnComp_alnRangeToSeqRange(struct stMalnComp *comp, int alnStart, int alnEnd, int *start, int *end) {
    struct stMalnCompCursor cc = stMalnCompCursor_make(comp);
    stMalnCompCursor_setAlignCol(&cc, alnStart);

    // find start
    if (!stMalnCompCursor_isAligned(&cc)) {
        if (!stMalnCompCursor_nextPos(&cc)) {
            *start = *end = 0;
            return false;
        }
    }
    *start = cc.pos;
    
    // find end
    stMalnCompCursor_setAlignCol(&cc,  alnEnd-1);
    *end = cc.pos;
    return true;
}


/* convert a sequence range to an alignment range, set range to 0-0 and return
 * false if none aligned */
bool stMalnComp_seqRangeToAlnRange(struct stMalnComp *comp, int start, int end, int *alnStart, int *alnEnd) {
    struct stMalnCompCursor cc = stMalnCompCursor_make(comp);
    stMalnCompCursor_setSeqPos(&cc, start);

    // find start
    if (!stMalnCompCursor_isAligned(&cc)) {
        if (!stMalnCompCursor_nextPos(&cc)) {
            *alnStart = *alnEnd = 0;
            return false;
        }
    }
    *alnStart = cc.alnIdx;
    
    // find end
    stMalnCompCursor_setSeqPos(&cc,  end-1);
    *alnEnd = cc.alnIdx+1;
    return true;
}

/* pad component to the specified alignment width */
void stMalnComp_pad(struct stMalnComp *comp, int width) {
    assert(stMalnComp_getWidth(comp) <= width);
    dyStringAppendMultiC(comp->alnStr, '-', width-stMalnComp_getWidth(comp));
}

/* append the specified number of characters from the string, adjusting component
 * coordinates */
void stMalnComp_append(struct stMalnComp *comp, char *src, int len) {
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
void stMalnComp_appendCompAln(struct stMalnComp *comp, struct stMalnComp *srcComp, int alnStart, int alnEnd) {
    stMalnComp_append(comp, stMalnComp_getAln(srcComp)+alnStart, alnEnd-alnStart);
}
