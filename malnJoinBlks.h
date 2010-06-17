#ifndef malnJoinBlks_h
#define malnJoinBlks_h
/* join two blocks using the specified reference components. */
struct malnBlk *malnJoinBlks(struct malnBlk *blk1, struct malnComp *refComp1, struct malnBlk *blk2, struct malnComp *refComp2);
#endif
