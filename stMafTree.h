#ifndef stMafTree_h
#define stMafTree_h
#include "sonLibTypes.h"
#include "mafJoinTypes.h"
struct mafAli;

/* Construct a MafTree object from an mafComp block if it contains a tree.  If
 * it doesn't contain a tree and is pairwise maf, create a tree with the root
 * the last component of the maf.  defaultBranchLength is used to assign
 * branch lengths when inferring trees from pair-wise MAFs. */
stMafTree *stMafTree_constructFromMaf(struct mafAli *ali, double defaultBranchLength);

/* clone a stMafTree */
stMafTree *stMafTree_clone(stMafTree *srcMTree);

/* destroy a MafTree object */
void stMafTree_destruct(stMafTree *mTree);

/* validate tree order and nodes matches mafAli */
void stMafTree_validateWithMaf(stMafTree *mTree, struct mafAli *ali);

/* compare function for coordinates by tree order  */
int stMafTree_treeOrderCmp(stMafTree *mTree, const char *orgSeq1, int chromStart1, int chromEnd1, const char *orgSeq2, int chromStart2, int chromEnd2);

/* join two trees at nodes specified by components */
stMafTree *stMafTree_join(stMafTree *mTree1, const char *orgSeq1, int chromStart1, int chromEnd1, stMafTree *mTree2, const char *orgSeq2, int chromStart2, int chromEnd2);

/* format tree as a newick string */
char *stMafTree_format(stMafTree *mTree);
#endif
