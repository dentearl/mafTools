#ifndef malnDeleteBlks_h
#define malnDeleteBlks_h
#include <stdbool.h>
struct malnBlk;

/* object for tracking blocks that are to be deleted */
struct malnDeleteBlks;

/* constructor */
struct malnDeleteBlks *malnDeleteBlks_construct(void);

/* delete blocks, removing from the set, and free the object */
void malnDeleteBlks_destruct(struct malnDeleteBlks *delBlks);

/* add a component's block to the delete table */
void malnDeleteBlks_flag(struct malnDeleteBlks *delBlks, struct malnBlk *blk);

/* is a block in the delete table? */
bool malnDeleteBlks_contains(struct malnDeleteBlks *delBlks, struct malnBlk *blk);

#endif
