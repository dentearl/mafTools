#ifndef malnComp_h
#define malnComp_h
#include <stdbool.h>
#include "genome.h"
#include "common.h"
#include "dystring.h"
#include "mafTree.h"

/* 
 * component of a multiple alignment block.
 */
struct malnComp {
    struct malnComp *next;
    struct malnBlk *blk;   // block we are associated with
    struct Seq *seq;       // genome sequence
    char strand;
    int start;             // start/end in strand coordinates
    int end;
    int chromStart;        // start/end in chrom coordinates
    int chromEnd;
    struct dyString *alnStr;
    struct mafTreeNodeCompLink *ncLink;
};

/* get a width of component */
static inline int malnComp_getWidth(struct malnComp *comp) {
    return comp->alnStr->stringSize;
}

/* get number of aligned bases */
static inline int malnComp_getAligned(struct malnComp *comp) {
    return (comp->chromEnd - comp->chromStart);
}

/* get a column */
static inline char malnComp_getCol(struct malnComp *comp, int alnIdx) {
    assert((0 <= alnIdx) && (alnIdx < comp->alnStr->stringSize));
    return comp->alnStr->string[alnIdx];
}

/* set a column */
static inline void malnComp_setCol(struct malnComp *comp, int alnIdx, char base) {
    assert((0 <= alnIdx) && (alnIdx < comp->alnStr->stringSize));
    comp->alnStr->string[alnIdx] = base;
}

/* get entire alignment string */
static inline char *malnComp_getAln(struct malnComp *comp) {
    return comp->alnStr->string;
}

/* is a position aligned? */
static inline bool malnComp_isAligned(struct malnComp *comp, int alnIdx) {
    return isBase(malnComp_getCol(comp, alnIdx));
}

/* component constructor */
struct malnComp *malnComp_construct(struct Seq *seq, char strand, int start, int end, char *alnStr);

/* component constructor to clone another component */
struct malnComp *malnComp_constructClone(struct malnComp *srcComp);

/* destructor */
void malnComp_destruct(struct malnComp *comp);

/* get tree location for component */
enum mafTreeLoc malnComp_getLoc(struct malnComp *comp);

/* component reverse complement */
struct malnComp *malnComp_reverseComplement(struct malnComp *comp);

/* convert a chrom range to a strand range for this component */
void malnComp_chromRangeToStrandRange(struct malnComp *comp, int chromStart, int chromEnd, int *start, int *end);

/* convert an alignment range to a sequence range, set range to 0-0 and return
 * false if none aligned */
bool malnComp_alnRangeToSeqRange(struct malnComp *comp, int alnStart, int alnEnd, int *start, int *end) ;

/* convert a sequence range to an alignment range, set range to 0-0 and return
 * false if none aligned */
bool malnComp_seqRangeToAlnRange(struct malnComp *comp, int start, int end, int *alnStart, int *alnEnd);

/* convert a sequence chrom range to an alignment range, set range to 0-0 and
 * return false if none aligned */
bool malnComp_seqChromRangeToAlnRange(struct malnComp *comp, int chromStart, int chromEnd, int *alnStart, int *alnEnd);

/* pad component to the specified alignment width */
void malnComp_pad(struct malnComp *comp, int width);

/* append the specified number of characters from the string, adjusting component
 * coordinates */
void malnComp_append(struct malnComp *comp, char *src, int len);

/* Append the specified alignment region of a component to this component.
 * The component must extend the range of this component. */
void malnComp_appendCompAln(struct malnComp *comp, struct malnComp *srcComp, int alnStart, int alnEnd);

/* get number of bases in component, not including inserts */
static inline bool malnComp_numBases(struct malnComp *comp) {
    return (comp->chromEnd - comp->chromStart);
}

/* compare two components to see if they overlap */
static inline bool malnComp_overlap(struct malnComp *comp, struct malnComp *comp2) {
    return (comp->seq == comp2->seq) && (comp->chromStart < comp2->chromEnd) && (comp->chromEnd > comp2->chromStart);
}

/* compare two components to see if they overlap or are adjacent */
static inline bool malnComp_overlapAdjacent(struct malnComp *comp, struct malnComp *comp2) {
    return (comp->seq == comp2->seq) && (comp->chromStart <= comp2->chromEnd) && (comp->chromEnd >= comp2->chromStart);
}

/* compare two components to see if they overlap or are adjacent on the same strand */
static inline bool malnComp_overlapAdjacentStrand(struct malnComp *comp, struct malnComp *comp2) {
    return (comp->seq == comp2->seq) && (comp->strand == comp2->strand) && (comp->start <= comp2->end) && (comp->end >= comp2->start);
}

/* compare to a range to see if they overlap */
static inline bool malnComp_overlapRange(struct malnComp *comp, struct Seq *seq, int chromStart, int chromEnd) {
    return (comp->seq == seq) && (comp->chromStart < chromEnd) && (comp->chromEnd > chromStart);
}

/* compare to a range to see if they overlap or are adjacent */
static inline bool malnComp_overlapAdjacentRange(struct malnComp *comp, struct Seq *seq, int chromStart, int chromEnd) {
    return (comp->seq == seq) && (comp->chromStart <= chromEnd) && (comp->chromEnd >= chromStart);
}

/* can two components be joined? */
bool malnComp_canJoin(struct malnComp *comp1, struct malnComp *comp2);

/* assert that the component is set-consistent */
void malnComp_assert(struct malnComp *comp);

/* compare two components for deterministic sorting */
int malnComp_cmp(struct malnComp *comp1, struct malnComp *comp2);

/* construct an component from a subrange of this component. Return NULL if
 * the subrange does not contain any aligned bases. */
struct malnComp *malnComp_constructSubrange(struct malnComp *comp, int alnStart, int alnEnd);

/* print base information describing a comp, newline not included */
void malnComp_prInfo(struct malnComp *comp, FILE *fh);

/* print a component for debugging purposes */
void malnComp_dumpv(struct malnComp *comp, FILE *fh, const char *label, va_list args);

/* print a component for debugging purposes */
void malnComp_dump(struct malnComp *comp, FILE *fh, const char *label, ...)
#if defined(__GNUC__)
__attribute__((format(printf, 3, 4)))
#endif
;

#endif
