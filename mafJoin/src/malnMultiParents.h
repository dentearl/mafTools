#ifndef malnMultiParents_h
#define malnMultiParents_h
#include <stdbool.h>
#include <stdio.h>
struct malnSet;
struct malnBlkSet;

/* check for regions with multiple parents, deleting blocks if if requested.
 * otherwise aborting. */
void malnMultiParents_check(struct malnSet *malnSet);
#endif
