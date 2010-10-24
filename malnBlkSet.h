#ifndef malnBlkSet_h
#define malnBlkSet_h
#include <stdbool.h>
struct malnBlk;

/* object to track blocks by unique id.  Linux x86_64 doesn't consistently
 * assign addresses for each run of a process, so we must have a unique
 * object id in each block */
struct malnBlkSet;

/* constructor */
struct malnBlkSet *malnBlkSet_construct(void);

/* make a clone of a block set, excluding ones marked as deleted */
struct malnBlkSet *malnBlkSet_constructClone(struct malnBlkSet *srcBlks);

/* destructor*/
void malnBlkSet_destruct(struct malnBlkSet *blks);

/* Add block to the map */
void malnBlkSet_add(struct malnBlkSet *blks, struct malnBlk *blk);

/* remove a block from the map */
void malnBlkSet_remove(struct malnBlkSet *blks, struct malnBlk *blk);

/* is a block in the map */
bool malnBlkSet_contains(struct malnBlkSet *blks, struct malnBlk *blk);

/* remove the first block from the map, which will include deleted blocks */
struct malnBlk *malnBlkSet_pop(struct malnBlkSet *blks);

/* get an iterator over non-deleted blocked in the map */
struct malnBlkSetIterator *malnBlkSet_getIterator(struct malnBlkSet *blks);

/* get the next block, NULL if no more */
struct malnBlk *malnBlkSetIterator_getNext(struct malnBlkSetIterator *iter);

/* destructor */
void malnBlkSetIterator_destruct(struct malnBlkSetIterator *iter);
#endif
