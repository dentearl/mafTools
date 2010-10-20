#ifndef malnJoin_h
#define malnJoin_h
struct Genome;
struct malnSet;

/* join two sets, generating a third */
struct malnSet *malnJoinSets(struct Genome *guideGenome, struct malnSet *malnSet1, struct malnSet *malnSet2, int maxBlkWidth);

#endif
