#include "malnBlkCursor.h"
#include "malnBlk.h"
#include "common.h"


/* construct a new cursor, if refComp is not null, then force it to be first */
struct malnBlkCursor *malnBlkCursor_construct(struct malnBlk *blk, struct malnComp *refComp) {
    // allocate as one memory block
    int numComps = slCount(blk->comps);
    struct malnBlkCursor *bc = needMem(sizeof(struct malnBlkCursor) + numComps*sizeof(struct malnCompCursor));
    bc->blk = blk;
    bc->rows = (struct malnCompCursor*)(((int8_t*)bc) + sizeof(struct malnBlkCursor));
    bc->numRows = numComps;
    int iRow = 0;

    // optionally force refComp first
    if (refComp != NULL) {
        bc->alnWidth = malnComp_getWidth(refComp);
        bc->rows[iRow++] = malnCompCursor_make(refComp);
    }
    
    for (struct malnComp *comp = blk->comps; comp != NULL; comp = comp->next) {
        if (comp != refComp) {
            if (iRow == 0) {
                bc->alnWidth = malnComp_getWidth(comp);
            }
            bc->rows[iRow++] = malnCompCursor_make(comp);
            assert(bc->alnWidth == malnComp_getWidth(comp));
        }
    }
    return bc;
}

/* destructor */
void malnBlkCursor_destruct(struct malnBlkCursor *bc) {
    freeMem(bc);
}
