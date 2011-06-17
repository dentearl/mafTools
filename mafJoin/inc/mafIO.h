#ifndef mafIO_h
#define mafIO_h
#include <stdbool.h>
#include <stdio.h>
struct Genome;
struct Genomes;
struct mafAli;

/* convert a mafAli to an malnBlk */
struct malnBlk *mafIO_malnBlkRead(struct Genomes *genomes, struct mafAli *mafAli, double defaultBranchLength, struct Genome *treelessRootGenome, bool requireTree);

/* Construct a malnSet from a MAF file. defaultBranchLength is used to
 * assign branch lengths when inferring trees from the MAF. */
struct malnSet *mafIO_malnSetRead(struct Genomes *genomes, char *mafFileName, int maxInputBlkWidth, double defaultBranchLength, struct Genome *treelessRootGenome);

/* write a block to a MAF */
void mafIO_malnBlkWrite(struct malnBlk *blk, FILE *mafFh);

/* write a malnSet to a MAF file  */
void mafIO_malnSetWrite(struct malnSet *malnSet, char *mafFileName);

#endif
