#ifndef malnCompCompMap_h
#define malnCompCompMap_h
struct malnComp;
#include <stdio.h>

/* map of component to another component. */
struct malnCompCompMap;

/* constructor */
struct malnCompCompMap *malnCompCompMap_construct(void);

/* destructor */
void malnCompCompMap_destruct(struct malnCompCompMap *mccm);

/* add a mapping */
void malnCompCompMap_add(struct malnCompCompMap *mccm, struct malnComp *srcComp, struct malnComp *destComp);

/* lookup a mapping, error if not found */
struct malnComp *malnCompCompMap_get(struct malnCompCompMap *mccm, struct malnComp *srcComp);

/* lookup a mapping, NULL in not found */
struct malnComp *malnCompCompMap_find(struct malnCompCompMap *mccm, struct malnComp *srcComp);

/* dump for debugging purposes */
void malnCompCompMap_dump(struct malnCompCompMap *mccm, FILE *fh, char *label);

#endif
