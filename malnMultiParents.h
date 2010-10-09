#ifndef malnMultiParents_h
#define malnMultiParents_h
#include <stdbool.h>
struct malnSet;

/* check for regions with multiple parents, deleting the blocks if requested,
 * otherwise aborting */
void malnMultiParents_check(struct malnSet *malnSet, bool discardTwoParents);
#endif
