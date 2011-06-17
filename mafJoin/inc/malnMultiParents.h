#ifndef malnMultiParents_h
#define malnMultiParents_h
#include <stdbool.h>
#include <stdio.h>
struct malnSet;
struct malnBlkSet;

/* open file for reporting problems, writing header */
FILE *malnMultiParents_openResolveDropLog(char *multiParentDroppedFile);

/* Resolve conflicts where blocks imply sequences with multiple parents,
 * editing blocks */
int malnMultiParents_resolve(struct malnSet *malnSet, FILE *dropLogFh);

/* check for regions with multiple parents, deleting blocks if if requested.
 * otherwise aborting. */
void malnMultiParents_check(struct malnSet *malnSet);
#endif
