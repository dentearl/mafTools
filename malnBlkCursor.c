#include "malnBlkCursor.h"
#include "malnBlk.h"
#include "common.h"

/* build cursor with all components */
static void constructWithAllComps(struct malnBlkCursor *bc, struct malnComp *refComp) {
    // optionally force refComp first
    if (refComp != NULL) {
        bc->alnWidth = malnComp_getWidth(refComp);
        bc->rows[bc->numRows++] = malnCompCursor_make(refComp);
    }
    
    for (struct malnComp *comp = bc->blk->comps; comp != NULL; comp = comp->next) {
        if (comp != refComp) {
            if (bc->numRows == 0) {
                bc->alnWidth = malnComp_getWidth(comp);
            }
            assert(bc->alnWidth == malnComp_getWidth(comp));
            bc->rows[bc->numRows++] = malnCompCursor_make(comp);
        }
    }
}

/* build cursor with a subset of components */
static void constructWithSubsetComps(struct malnBlkCursor *bc, struct malnComp **subsetComps) {
    while (subsetComps[bc->numRows] != NULL) {
        assert(subsetComps[bc->numRows]->blk == bc->blk);
        if (bc->numRows == 0) {
            bc->alnWidth = malnComp_getWidth(subsetComps[bc->numRows]);
        }
        assert(bc->alnWidth == malnComp_getWidth(subsetComps[bc->numRows]));
        bc->rows[bc->numRows] = malnCompCursor_make(subsetComps[bc->numRows]);
        bc->numRows++;
    }
}

/* construct a new cursor, if refComp is not null, then force it to be first. If subsetComps
 * is not NULL, it is a NULL terminated list of a subset of components to include
 * in the order specified. */
struct malnBlkCursor *malnBlkCursor_construct(struct malnBlk *blk, struct malnComp *refComp, struct malnComp **subsetComps) {
    assert((refComp == NULL) || (subsetComps == NULL) || (refComp == subsetComps[0]));

    // allocate as one memory block, assume all components are included
    struct malnBlkCursor *bc = needMem(sizeof(struct malnBlkCursor) + slCount(blk->comps)*sizeof(struct malnCompCursor));
    bc->blk = blk;
    bc->rows = (struct malnCompCursor*)(((int8_t*)bc) + sizeof(struct malnBlkCursor));
    if (subsetComps == NULL) {
        constructWithAllComps(bc, refComp);
    } else {
        constructWithSubsetComps(bc, subsetComps);
    }
    return bc;
}

/* destructor */
void malnBlkCursor_destruct(struct malnBlkCursor *bc) {
    freeMem(bc);
}
