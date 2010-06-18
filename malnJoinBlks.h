#ifndef malnJoinBlks_h
#define malnJoinBlks_h
struct malnComp;

/* join two blocks using their specified reference components. */
struct malnBlk *malnJoinBlks(struct malnComp *refComp1, struct malnComp *refComp2);
#endif
