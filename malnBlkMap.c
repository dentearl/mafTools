#include "malnBlkMap.h"
#include "sonLibSortedSet.h"
#include "malnBlk.h"
#include "common.h"

struct malnBlkMap {
    stSortedSet *map;
};

struct malnBlkMapIterator {
    stSortedSetIterator *iter;
};

/* compare function for block objId */
static int objIdCmp(const void *vblk1, const void *vblk2) {
    struct malnBlk *blk1 = (struct malnBlk *)vblk1;
    struct malnBlk *blk2 = (struct malnBlk *)vblk2;
    return blk1->objId - blk2->objId;
}


struct malnBlkMap *malnBlkMap_construct(void) {
    struct malnBlkMap *blks;
    AllocVar(blks);
    blks->map = stSortedSet_construct3(objIdCmp, NULL);
    return blks;
}

void malnBlkMap_destruct(struct malnBlkMap *blks) {
    if (blks != NULL) {
        stSortedSet_destruct(blks->map);
        freeMem(blks);
    }
}

void malnBlkMap_add(struct malnBlkMap *blks, struct malnBlk *blk) {
    stSortedSet_insert(blks->map, blk);
}

void malnBlkMap_remove(struct malnBlkMap *blks, struct malnBlk *blk) {
    stSortedSet_remove(blks->map, blk);
}

bool malnBlkMap_contains(struct malnBlkMap *blks, struct malnBlk *blk) {
    return stSortedSet_search(blks->map, blk) != NULL; 
}

void malnBlkMap_deleteAll(struct malnBlkMap *blks) {
    struct malnBlk *blk;
    while ((blk = stSortedSet_getFirst(blks->map)) != NULL) {
        stSortedSet_remove(blks->map, blk);
        malnBlk_destruct(blk);
    }
}

struct malnBlkMapIterator *malnBlkMap_getIterator(struct malnBlkMap *blks) {
    struct malnBlkMapIterator *iter;
    AllocVar(iter);
    iter->iter = stSortedSet_getIterator(blks->map);
    return iter;
}

struct malnBlk *malnBlkMapIterator_getNext(struct malnBlkMapIterator *iter) {
    return stSortedSet_getNext(iter->iter);
}

void malnBlkMapIterator_destruct(struct malnBlkMapIterator *iter) {
    stSortedSet_destructIterator(iter->iter);
    freeMem(iter);
}
