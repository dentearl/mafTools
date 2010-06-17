#ifndef stMalnComp_h
#define stMalnComp_h
#include <stdbool.h>
#include "genome.h"
#include "common.h"
#include "dystring.h"

/* location of a component in the tree */
enum stMalnCompTreeLoc {
    stMalnCompTreeUnknown,
    stMalnCompTreeRoot,
    stMalnCompTreeInternal,
    stMalnCompTreeLeaf
};


/* 
 * component of a multiple alignment block.
 */
struct stMalnComp {
    struct stMalnComp *next;
    struct Seq *seq;       // genome sequence
    char strand;
    int start;             // start/end in strand coordinates
    int end;
    int chromStart;        // start/end in chrom coordinates
    int chromEnd;
    struct dyString *alnStr;
    bool isReference;      // is this the reference genome?
    enum stMalnCompTreeLoc treeLoc;  // location in tree.
};

/* get a width of component */
static inline int stMalnComp_getWidth(struct stMalnComp *comp) {
    return comp->alnStr->stringSize;
}

/* get a column */
static inline char stMalnComp_getCol(struct stMalnComp *comp, int alnIdx) {
    assert((0 <= alnIdx) && (alnIdx < comp->alnStr->stringSize));
    return comp->alnStr->string[alnIdx];
}

/* get entire alignment string */
static inline char *stMalnComp_getAln(struct stMalnComp *comp) {
    return comp->alnStr->string;
}

/* is a position aligned? */
static inline bool stMalnComp_isAligned(struct stMalnComp *comp, int alnIdx) {
    return isBase(stMalnComp_getCol(comp, alnIdx));
}

/* component constructor */
struct stMalnComp *stMalnComp_construct(struct Seq *seq, char strand, int start, int end, char *alnStr);

/* component constructor to clone another component */
struct stMalnComp *stMalnComp_constructClone(struct stMalnComp *srcComp);

/* destructor */
void stMalnComp_destruct(struct stMalnComp *comp);

/* component reverse complement */
struct stMalnComp *stMalnComp_reverseComplement(struct stMalnComp *comp);

/* convert an alignment range to a sequence range, set range to 0-0 and return
 * false if none aligned */
bool stMalnComp_alnRangeToSeqRange(struct stMalnComp *comp, int alnStart, int alnEnd, int *start, int *end) ;

/* convert a sequence range to an alignment range, set range to 0-0 and return
 * false if none aligned */
bool stMalnComp_seqRangeToAlnRange(struct stMalnComp *comp, int start, int end, int *alnStart, int *alnEnd);

/* pad component to the specified alignment width */
void stMalnComp_pad(struct stMalnComp *comp, int width);

/* append the specified number of characters from the string, adjusting component
 * coordinates */
void stMalnComp_append(struct stMalnComp *comp, char *src, int len);

/* Append the specified alignment region of a component to this component.
 * The component must extend the range of this component. */
void stMalnComp_appendCompAln(struct stMalnComp *comp, struct stMalnComp *srcComp, int alnStart, int alnEnd);

/* compare to component to see if they overlap */
static inline bool stMalnComp_overlap(struct stMalnComp *comp, struct stMalnComp *comp2) {
    return (comp->seq == comp2->seq) && (comp->chromStart < comp2->chromEnd) && (comp->chromEnd > comp2->chromStart);
}

#endif
