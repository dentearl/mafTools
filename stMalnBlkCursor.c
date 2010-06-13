#include "stMalnBlkCursor.h"
#include "stMalnBlk.h"
#include "common.h"


/* construct a new cursor, if refComp is not null, then force it to be first */
struct stMalnBlkCursor *stMalnBlkCursor_construct(struct stMalnBlk *blk, struct stMalnComp *refComp) {
    // allocate as one memory block
    int numComps = slCount(blk->comps);
    struct stMalnBlkCursor *bc = needMem(sizeof(struct stMalnBlkCursor) + numComps*sizeof(struct stMalnCompCursor));
    bc->blk = blk;
    bc->rows = (struct stMalnCompCursor*)(((int8_t*)bc) + sizeof(struct stMalnBlkCursor));
    bc->numRows = numComps;
    int iRow = 0;

    // optionally force refComp first
    if (refComp != NULL) {
        bc->alnWidth = stMalnComp_getWidth(refComp);
        bc->rows[iRow++] = stMalnCompCursor_make(refComp);
    }
    
    for (struct stMalnComp *comp = blk->comps; comp != NULL; comp = comp->next) {
        if (comp != refComp) {
            if (iRow == 0) {
                bc->alnWidth = stMalnComp_getWidth(comp);
            }
            bc->rows[iRow++] = stMalnCompCursor_make(comp);
            assert(bc->alnWidth == stMalnComp_getWidth(comp));
        }
    }
    return bc;
}

/* destructor */
void stMalnBlkCursor_destruct(struct stMalnBlkCursor *bc) {
    freeMem(bc);
}
