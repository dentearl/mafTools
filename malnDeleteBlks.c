#include "malnDeleteBlks.h"
#include "sonLibHash.h"
#include "malnBlk.h"
#include "common.h"

struct malnDeleteBlks {
    stHash *blks;
};

struct malnDeleteBlks *malnDeleteBlks_construct(void) {
    struct malnDeleteBlks *delBlks;
    AllocVar(delBlks);
    delBlks->blks = stHash_construct2(NULL, (void (*)(void *))malnBlk_destruct);
    return delBlks;
}

void malnDeleteBlks_destruct(struct malnDeleteBlks *delBlks) {
    if (delBlks != NULL) {
        stHash_destruct(delBlks->blks);
        freeMem(delBlks);
    }
}

/* add a component's block to the delete table */
void malnDeleteBlks_flag(struct malnDeleteBlks *delBlks, struct malnBlk *blk) {
    stHash_insert(delBlks->blks, blk, blk);
}

/* is a block in the delete table? */
bool malnDeleteBlks_contains(struct malnDeleteBlks *delBlks, struct malnBlk *blk) {
    return stHash_search(delBlks->blks, blk) != NULL; 
}
