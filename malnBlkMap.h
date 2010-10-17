#ifndef malnBlkMap_h
#define malnBlkMap_h
#include <stdbool.h>
struct malnBlk;

/* object to track blocks by unique id.  Linux x86_64 doesn't consistently
 * assign addresses for each run of a process, so we must have a unique
 * object id in each block */
struct malnBlkMap;

/* constructor */
struct malnBlkMap *malnBlkMap_construct(void);

/* destructor*/
void malnBlkMap_destruct(struct malnBlkMap *blks);

/* Add block to the map */
void malnBlkMap_add(struct malnBlkMap *blks, struct malnBlk *blk);

/* remove a block from the map */
void malnBlkMap_remove(struct malnBlkMap *blks, struct malnBlk *blk);

/* is a block in the map */
bool malnBlkMap_contains(struct malnBlkMap *blks, struct malnBlk *blk);

/* remove the first block from the map */
struct malnBlk *malnBlkMap_pop(struct malnBlkMap *blks);

/* get an iterator over the items in the map */
struct malnBlkMapIterator *malnBlkMap_getIterator(struct malnBlkMap *blks);

/* get the next block, NULL if no more */
struct malnBlk *malnBlkMapIterator_getNext(struct malnBlkMapIterator *iter);

/* destructor */
void malnBlkMapIterator_destruct(struct malnBlkMapIterator *iter);
#endif
