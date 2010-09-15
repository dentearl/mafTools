#ifndef malnComp_h
#define malnComp_h
#include <stdbool.h>
#include "genome.h"
#include "common.h"
#include "dystring.h"

/* location of a component in the tree, can also be a bit-set for
 * selection. */
enum malnCompTreeLoc {
    malnCompTreeUnknown    = 0x0,
    malnCompTreeRoot       = 0x1,
    malnCompTreeInternal   = 0x2,
    malnCompTreeLeaf       = 0x4
};


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
    enum malnCompTreeLoc treeLoc;  // location in tree.
};

/* get a width of component */
static inline int malnComp_getWidth(struct malnComp *comp) {
    return comp->alnStr->stringSize;
}

/* get a column */
static inline char malnComp_getCol(struct malnComp *comp, int alnIdx) {
    assert((0 <= alnIdx) && (alnIdx < comp->alnStr->stringSize));
    return comp->alnStr->string[alnIdx];
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

/* component reverse complement */
struct malnComp *malnComp_reverseComplement(struct malnComp *comp);

/* convert an alignment range to a sequence range, set range to 0-0 and return
 * false if none aligned */
bool malnComp_alnRangeToSeqRange(struct malnComp *comp, int alnStart, int alnEnd, int *start, int *end) ;

/* convert a sequence range to an alignment range, set range to 0-0 and return
 * false if none aligned */
bool malnComp_seqRangeToAlnRange(struct malnComp *comp, int start, int end, int *alnStart, int *alnEnd);

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

/* compare to component to see if they overlap */
static inline bool malnComp_overlap(struct malnComp *comp, struct malnComp *comp2) {
    return (comp->seq == comp2->seq) && (comp->chromStart < comp2->chromEnd) && (comp->chromEnd > comp2->chromStart);
}

/* compare to a range to see if they overlap */
static inline bool malnComp_overlapRange(struct malnComp *comp, struct Seq *seq, int chromStart, int chromEnd) {
    return (comp->seq == seq) && (comp->chromStart < chromEnd) && (comp->chromEnd > chromStart);
}

/* can a join be made at this component? */
static inline bool malnComp_joinable(struct malnComp *comp) {
    return (comp->treeLoc & (malnCompTreeRoot|malnCompTreeLeaf)) != 0;
}

/* can two components be joined? */
static inline bool malnComp_canJoin(struct malnComp *comp1, struct malnComp *comp2) {
    return malnComp_overlap(comp1, comp2)
        && (((comp1->treeLoc == malnCompTreeRoot) && (comp2->treeLoc == malnCompTreeRoot))
            || ((comp1->treeLoc == malnCompTreeRoot) && (comp2->treeLoc == malnCompTreeLeaf))
            || ((comp1->treeLoc == malnCompTreeLeaf) && (comp2->treeLoc == malnCompTreeRoot)));
}

/* assert that the component is set-consistent */
void malnComp_assert(struct malnComp *comp);

/* print a component for debugging purposes */
void malnComp_dump(struct malnComp *comp, const char *label, FILE *fh);

#endif
