#include "malnBlkSet.h"
#include "sonLibSortedSet.h"
#include "malnBlk.h"
#include "common.h"

struct malnBlkSet {
    stSortedSet *map;
};

struct malnBlkSetIterator {
    stSortedSetIterator *iter;
};

/* compare function for block objId */
static int objIdCmp(const void *vblk1, const void *vblk2) {
    struct malnBlk *blk1 = (struct malnBlk *)vblk1;
    struct malnBlk *blk2 = (struct malnBlk *)vblk2;
    return blk1->objId - blk2->objId;
}


struct malnBlkSet *malnBlkSet_construct(void) {
    struct malnBlkSet *blks;
    AllocVar(blks);
    blks->map = stSortedSet_construct3(objIdCmp, NULL);
    return blks;
}

void malnBlkSet_destruct(struct malnBlkSet *blks) {
    if (blks != NULL) {
        stSortedSet_destruct(blks->map);
        freeMem(blks);
    }
}

void malnBlkSet_add(struct malnBlkSet *blks, struct malnBlk *blk) {
    stSortedSet_insert(blks->map, blk);
}

void malnBlkSet_remove(struct malnBlkSet *blks, struct malnBlk *blk) {
    stSortedSet_remove(blks->map, blk);
}

bool malnBlkSet_contains(struct malnBlkSet *blks, struct malnBlk *blk) {
    return stSortedSet_search(blks->map, blk) != NULL; 
}

struct malnBlk *malnBlkSet_pop(struct malnBlkSet *blks) {
    struct malnBlk *blk = stSortedSet_getFirst(blks->map);
    if (blk != NULL) {
        stSortedSet_remove(blks->map, blk);
    }
    return blk;
}

struct malnBlkSetIterator *malnBlkSet_getIterator(struct malnBlkSet *blks) {
    struct malnBlkSetIterator *iter;
    AllocVar(iter);
    iter->iter = stSortedSet_getIterator(blks->map);
    return iter;
}

struct malnBlk *malnBlkSetIterator_getNext(struct malnBlkSetIterator *iter) {
    return stSortedSet_getNext(iter->iter);
}

void malnBlkSetIterator_destruct(struct malnBlkSetIterator *iter) {
    stSortedSet_destructIterator(iter->iter);
    freeMem(iter);
}
