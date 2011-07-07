#ifndef malnComp_h
#define malnComp_h
#include <stdbool.h>
#include "genome.h"
#include "common.h"
#include "mafTree.h"
struct malnCompCursor;

/*
 * component segment: a contiguous range of bases in a component of an alignment block.
 */
struct malnCompSeg {
    struct malnCompSeg *next;  // link to next segment, in ascending order
    int start;     // start of region in sequence strand coordinates
    int alnStart;  // start of region in alignment
    int size;      // size of region
    char *bases;   // sequence bases
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
    int alnWidth;          // width of this sequence alignment
    struct malnCompSeg *segs;
    struct mafTreeNodeCompLink *ncLink;
};

/* clone a segment */
struct malnCompSeg *malnCompSeg_constructClone(struct malnCompSeg *srcSeg);

/* extend a segment */
void malnCompSeg_extend(struct malnCompSeg *seg, char *newBases, int newSize);

/* get sequence end position for a seg */
static inline int malnCompSeg_getEnd(struct malnCompSeg *seg) {
    return seg->start + seg->size;
}

/* get alignment end position for a seg */
static inline int malnCompSeg_getAlnEnd(struct malnCompSeg *seg) {
    return seg->alnStart + seg->size;
}

/* get a base from a seg */
static inline char malnCompSeg_getBase(struct malnCompSeg *seg, int alnIdx) {
    assert((seg->alnStart <= alnIdx) && (alnIdx < malnCompSeg_getAlnEnd(seg)));
    return seg->bases[alnIdx - seg->alnStart];
}

/* is the specified alignment index in the seg */
static inline bool malnCompSeg_alnIdxInSeg(struct malnCompSeg *seg, int alnIdx) {
    return (seg->alnStart <= alnIdx) && (alnIdx < malnCompSeg_getAlnEnd(seg));
}

/* is the specified position in the seg */
static inline bool malnCompSeg_containsPos(struct malnCompSeg *seg, int pos) {
    return (seg->start <= pos) && (pos < malnCompSeg_getEnd(seg));
}

/* Are two segments overlapping or adjacent in position space.  segments maybe
 * in different components. */
static inline bool malnCompSeg_overlapAdjacentStrand(struct malnCompSeg *seg, struct malnComp *comp,
                                                     struct malnCompSeg *seg2, struct malnComp *comp2) {
    return (comp->seq == comp2->seq) && (comp->strand == comp2->strand)
        && (seg->start <= malnCompSeg_getEnd(seg2))
        && (malnCompSeg_getEnd(seg) >= seg2->start);
}

/* get number of aligned bases */
static inline int malnComp_getNumAligned(struct malnComp *comp) {
    return (comp->chromEnd - comp->chromStart);
}

/* constructor */
struct malnComp *malnComp_construct(struct Seq *seq, char strand, int start, int end, int alnWidth);

/* constructor from alignment string  */
struct malnComp *malnComp_constructFromStr(struct Seq *seq, char strand, int start, int end, char *alnStr);

/* component constructor to clone another component */
struct malnComp *malnComp_constructClone(struct malnComp *srcComp);

/* destructor */
void malnComp_destruct(struct malnComp *comp);

/* Clear sequence to reduce memory associated with dying blocks */
void malnComp_freeSeqMem(struct malnComp *comp);

/* get tree location for component */
enum mafTreeLoc malnComp_getLoc(struct malnComp *comp);

/* component reverse complement */
struct malnComp *malnComp_reverseComplement(struct malnComp *srcComp);

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
static inline void malnComp_pad(struct malnComp *comp, int width) {
    assert(comp->alnWidth <= width);
    comp->alnWidth = width;
}

/* Append from the current position of the specified cursor up to the
 * given alignment column   */
void malnComp_appendFromCursor(struct malnComp *comp, struct malnCompCursor *srcCursor, int alnEnd);

/* get number of bases in component, not including inserts */
static inline bool malnComp_numAlignedBases(struct malnComp *comp) {
    return (comp->chromEnd - comp->chromStart);
}

/* compare two components to see if they overlap */
static inline bool malnComp_overlap(struct malnComp *comp, struct malnComp *comp2) {
    return (comp->seq == comp2->seq) && (comp->chromStart < comp2->chromEnd) && (comp->chromEnd > comp2->chromStart);
}

/* compare two components to see if they overlap or are adjacent on the same strand */
static inline bool malnComp_overlapStrand(struct malnComp *comp, struct malnComp *comp2) {
    return (comp->seq == comp2->seq) && (comp->strand == comp2->strand) && (comp->start < comp2->end) && (comp->end > comp2->start);
}

/* compare two components to see if they overlap or are adjacent */
static inline bool malnComp_overlapAdjacent(struct malnComp *comp, struct malnComp *comp2) {
    return (comp->seq == comp2->seq) && (comp->chromStart <= comp2->chromEnd) && (comp->chromEnd >= comp2->chromStart);
}

/* compare two components to see if they overlap or are adjacent on the same strand */
static inline bool malnComp_overlapAdjacentStrand(struct malnComp *comp, struct malnComp *comp2) {
    return (comp->seq == comp2->seq) && (comp->strand == comp2->strand) && (comp->start <= comp2->end) && (comp->end >= comp2->start);
}

/* compare two components to see if they are adjacent but not overlapped on the same strand */
static inline bool malnComp_adjacentStrand(struct malnComp *comp, struct malnComp *comp2) {
    return (comp->seq == comp2->seq) && (comp->strand == comp2->strand) && ((comp->start == comp2->end) || (comp->end == comp2->start));
}

/* compare to a range to see if they overlap */
static inline bool malnComp_overlapRange(struct malnComp *comp, struct Seq *seq, int chromStart, int chromEnd) {
    return (comp->seq == seq) && (comp->chromStart < chromEnd) && (comp->chromEnd > chromStart);
}

/* compare to a range to see if they overlap or are adjacent */
static inline bool malnComp_overlapAdjacentRange(struct malnComp *comp, struct Seq *seq, int chromStart, int chromEnd) {
    return (comp->seq == seq) && (comp->chromStart <= chromEnd) && (comp->chromEnd >= chromStart);
}

/* compare two components to see if they overlap or are adjacent when the second 
 * component is in the specified relative orientation (-1 == reverse complemented) */
bool malnComp_overlapAdjacentOrient(struct malnComp *comp, struct malnComp *comp2, int orient);

/* is a stand position in the component */
static inline bool malnComp_containsPos(struct malnComp *comp, int pos) {
    return (comp->start <= pos) && (pos < comp->end);
}

/* can two components be joined? */
bool malnComp_canJoin(struct malnComp *comp1, struct malnComp *comp2);

/* assert that the component is set-consistent. ncLinkCheck can be false
 * for sanity check before links are constructed. */
void malnComp_assert(struct malnComp *comp, bool ncLinkCheck);

/* compare two components for deterministic sorting */
int malnComp_cmp(struct malnComp *comp1, struct malnComp *comp2);

/* compare two components in chrom order */
int malnComp_chromCmp(struct malnComp *comp1, struct malnComp *comp2);

/* construct an component from a subrange of this component. Return NULL if
 * the subrange does not contain any aligned bases. */
struct malnComp *malnComp_constructSubrange(struct malnComp *srcComp, int alnStart, int alnEnd);

/* are there any aligned bases in a range */
bool malnComp_anyAlignedRange(struct malnComp *comp, int alnStart, int alnEnd);

/* get the last block in the component */
static inline struct malnCompSeg *malnComp_getLastSeg(struct malnComp *comp) {
    struct malnCompSeg *seg = comp->segs;
    while (seg->next != NULL) {
        seg = seg->next;
    }
    return seg;
}

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
