#ifndef malnCompCompMap_h
#define malnCompCompMap_h
struct malnComp;

/* map of component to another component. */
struct malnCompCompMap;

/* constructor */
struct malnCompCompMap *malnCompCompMap_construct(void);

/* destructor */
void malnCompCompMap_destruct(struct malnCompCompMap *mccm);

/* insert a mapping */
void malnCompCompMap_insert(struct malnCompCompMap *mccm, struct malnComp *srcComp, struct malnComp *destComp);

/* lookup a mapping, error if not found */
struct malnComp *malnCompCompMap_get(struct malnCompCompMap *mccm, struct malnComp *srcComp);

#endif
