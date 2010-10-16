#ifndef malnJoinDups_h
#define malnJoinDups_h
#include <stdbool.h>
struct malnSet;

/* Join duplication blocks in a set, which evolver outputs as separate
 * blocks.*/
void malnJoinDups_joinSetDups(struct malnSet *malnSet);
#endif
