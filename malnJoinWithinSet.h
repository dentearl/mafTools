#ifndef malnJoinWithinSet_h
#define malnJoinWithinSet_h
#include <stdbool.h>
struct malnSet;

/* Join duplication blocks in a set, which evolver outputs as separate
 * blocks.*/
void malnJoinWithinSet_joinDups(struct malnSet *malnSet);

/* Join adjacent and overlapping blocks in a set.  Stop joining at columns
 * where only the root component crosses that column and the root is adjacent
 * and not overlapping. */
void malnJoinWithinSet_joinOverlapAdjacent(struct malnSet *malnSet);

#endif
