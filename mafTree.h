#ifndef mafTree_h
#define mafTree_h
#include "sonLibTypes.h"
#include "mafJoinTypes.h"
struct mafAli;
struct malnComp;
struct malnBlk;
struct Genome;
struct Genomes;

/* location of a component in the tree */
enum mafTreeLoc {
    mafTreeRoot,
    mafTreeInternal,
    mafTreeLeaf
};

/* Construct a MafTree object from an mafComp block if it contains a tree.  If
 * it doesn't contain a tree, create one with the root the last component of
 * the MAF or treelessRootGenome if specified.  defaultBranchLength is used to
 * assign branch lengths when inferring trees from the MAFs. */
mafTree *mafTree_constructFromMaf(struct Genomes *genomes, struct mafAli *ali, double defaultBranchLength, struct Genome *treelessRootGenome);

/* clone a mafTree */
mafTree *mafTree_clone(mafTree *srcMTree);

/* destroy a MafTree object */
void mafTree_destruct(mafTree *mTree);

/* validate tree order and nodes matches mafAli */
void mafTree_validateWithMaf(mafTree *mTree, struct mafAli *ali);

/* compare function for component coordinates by tree order  */
int mafTree_treeOrderCmp(mafTree *mTree, const char *orgSeq1, int chromStart1, int chromEnd1, const char *orgSeq2, int chromStart2, int chromEnd2);

/* join two trees at nodes specified by components */
mafTree *mafTree_join(mafTree *mTree1, const char *orgSeq1, int chromStart1, int chromEnd1, mafTree *mTree2, const char *orgSeq2, int chromStart2, int chromEnd2);

/* get location type in tree */
enum mafTreeLoc mafTree_getLoc(mafTree *mTree, const char *orgSeq, int chromStart, int chromEnd);

/* format tree as a newick string */
char *mafTree_format(mafTree *mTree);
#endif
