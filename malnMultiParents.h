#ifndef malnMultiParents_h
#define malnMultiParents_h
#include <stdbool.h>
struct malnSet;
struct malnBlkSet;

/* Resolve conflicts where blocks imply sequences with multiple parents,
 * editing blocks */
void malnMultiParents_resolve(struct malnSet *malnSet);

/* check for regions with multiple parents, deleting blocks if if requested.
 * otherwise aborting. */
void malnMultiParents_check(struct malnSet *malnSet);
#endif
