#ifndef malnBlk_h
#define malnBlk_h
#include "mafJoinTypes.h"
#include <stdbool.h>
struct Seq;
struct Genome;

/* 
 * Multiple alignment block.
 */
struct stMalnBlk {
    struct stMalnBlk *next;
    int alnWidth;
    struct stMalnComp *comps;   // components
    stMafTree *mTree;           // tree associated with block
    bool done;                  // finished some kind of processing
};

/* constructor */
struct stMalnBlk *stMalnBlk_construct(stMafTree *mTree);

/* constructor clone */
struct stMalnBlk *stMalnBlk_constructClone(struct stMalnBlk *srcBlk);

/* destructor */
void stMalnBlk_destruct(struct stMalnBlk *blk);

/* set location and type attributes from tree */
void stMalnBlk_setLocAttr(struct stMalnBlk *blk, struct Genome *refGenome);

/* add a component */
void stMalnBlk_addComp(struct stMalnBlk *blk, struct stMalnComp *comp);

/* sort components by tree */
void stMalnBlk_sortComps(struct stMalnBlk *blk);

/* find a component by seq and start, NULL if not found  */
struct stMalnComp *stMalnBlk_findBySeqStart(struct stMalnBlk *blk, struct Seq *seq, int start);

/* find a component by seq and chrom range, NULL if not found  */
struct stMalnComp *stMalnBlk_findCompByChromRange(struct stMalnBlk *blk, struct Seq *seq, int chromStart, int chromEnd);

/* block reverse complement */
struct stMalnBlk *stMalnBlk_reverseComplement(struct stMalnBlk *blk);

/* pad all components so they are the same width as the overall alignment. */
void stMalnBlk_pad(struct stMalnBlk *blk);

/* assert that the block is set-consistent */
void stMalnBlk_assert(struct stMalnBlk *blk);

#endif
