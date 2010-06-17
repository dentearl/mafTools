#ifndef pairBranchLength_h
#define pairBranchLength_h
#include "sonLibTypes.h"

error: do not use this

/* get distance between two species in the tree as an estimate of length between
 * species. */
double pairBranchLength(ETree *initBranchLenTree, const char *org1, const char *org2);
#endif
