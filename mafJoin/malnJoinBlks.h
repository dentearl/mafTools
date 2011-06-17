#ifndef malnJoinBlks_h
#define malnJoinBlks_h
struct malnComp;

/* join two blocks using their specified guide components. */
struct malnBlk *malnJoinBlks(struct malnComp *guideComp1, struct malnComp *guideComp2);
#endif
