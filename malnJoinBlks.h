#ifndef malnJoinBlks_h
#define malnJoinBlks_h
struct malnComp;

/* join two blocks using their specified reference components.  Optionally return
 * resulting join component. */
struct malnBlk *malnJoinBlks(struct malnComp *refComp1, struct malnComp *refComp2, struct malnComp **joinedCompRet);
#endif
