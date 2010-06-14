#ifndef stMalnSet_h
#define stMalnSet_h
#include "sonLibTypes.h"
struct Genomes;
struct Genome;
struct stMalnSet;
struct stMalnBlk;
struct Seq;

/* construct an empty stMalnSet  */
struct stMalnSet *stMalnSet_construct(struct Genomes *genomes, struct Genome *refGenome);

/* Construct a stMalnSet from a MAF file. defaultBranchLength is used to
 * assign branch lengths when inferring trees from pair-wise MAFs. */
struct stMalnSet *stMalnSet_constructFromMaf(struct Genomes *genomes, struct Genome *refGenome, char *mafFileName, double defaultBranchLength);

/* get associated genomes object  */
struct Genomes *stMalnSet_getGenomes(struct stMalnSet *malnSet);

/* add a block to a malnSet */
void stMalnSet_addBlk(struct stMalnSet *malnSet, struct stMalnBlk *blk);

/* write a stMalnSet to a MAF file  */
void stMalnSet_writeMaf(struct stMalnSet *malnSet, char *mafFileName);

/* get list of blocks.  Don't modify this list. */
struct stMalnBlk *stMalnSet_getBlocks(struct stMalnSet *malnSet);

/* get list of slRefs to blks who's reference range overlaps the specified
 * range.  slFreeList of results when done. */
struct slRef *stMalnSet_getOverlapping(struct stMalnSet *malnSet, struct Seq *seq, int chromStart, int chromEnd);

#endif
