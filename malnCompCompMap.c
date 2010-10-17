#include "malnCompCompMap.h"
#include "common.h"
#include "sonLibHash.h"


struct malnCompCompMap {
    stHash *map;
};

/* constructor */
struct malnCompCompMap *malnCompCompMap_construct(void) {
    struct malnCompCompMap *mccm;
    AllocVar(mccm);
    mccm->map = stHash_construct();
    return mccm;
}

/* destructor */
void malnCompCompMap_destruct(struct malnCompCompMap *mccm) {
    if (mccm != NULL) {
        stHash_destruct(mccm->map);
        freeMem(mccm);
    }
}

/* add a mapping */
void malnCompCompMap_add(struct malnCompCompMap *mccm, struct malnComp *srcComp, struct malnComp *destComp) {
    stHash_insert(mccm->map, srcComp, destComp);
}

/* lookup a mapping, error if not found */
struct malnComp *malnCompCompMap_get(struct malnCompCompMap *mccm, struct malnComp *srcComp) {
    struct malnComp *destComp = stHash_search(mccm->map, srcComp);
    if (destComp == NULL) {
        errAbort("no mapping found for component");
    }
    return destComp;
}

/* lookup a mapping, NULL in not found */
struct malnComp *malnCompCompMap_find(struct malnCompCompMap *mccm, struct malnComp *srcComp) {
    return stHash_search(mccm->map, srcComp);
}
