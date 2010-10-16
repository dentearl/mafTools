#ifndef malnMultiParents_h
#define malnMultiParents_h
#include <stdbool.h>
struct malnSet;
struct malnBlkMap;

/* check for regions with multiple parents, deleting the blocks if requested,
 * otherwise aborting.  If deletingMap is not, skip entries in it.  */
void malnMultiParents_check(struct malnSet *malnSet, bool discardTwoParents, struct malnBlkMap *deletedBlks);
#endif
