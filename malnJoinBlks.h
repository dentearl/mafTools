#ifndef malnJoinBlks_h
#define malnJoinBlks_h
struct malnComp;

/* join two blocks using their specified guide components.  Optionally return
 * resulting join component. */
struct malnBlk *malnJoinBlks(struct malnComp *guideComp1, struct malnComp *guideComp2, struct malnComp **joinedCompRet);
#endif
