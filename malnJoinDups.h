#ifndef malnJoinDups_h
#define malnJoinDups_h
#include <stdbool.h>
struct malnSet;

/* Join duplication blocks in a set, which evolver outputs as separate
 * blocks.  Duplications will only be joined at the root */
void malnJoin_joinSetDups(struct malnSet *malnSet, bool discardTwoParents);
#endif
