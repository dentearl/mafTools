#ifndef stMalnJoin_h
#define stMalnJoin_h
struct Genome;
struct stMalnSet;

/* join two sets, generating a third */
struct stMalnSet *stMalnJoin_joinSets(struct Genome *refGenome, struct stMalnSet *malnSet1, struct stMalnSet *malnSet2);

#endif
