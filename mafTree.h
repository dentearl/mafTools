#ifndef mafTree_h
#define mafTree_h
#include "sonLibTypes.h"
#include "mafJoinTypes.h"
struct malnCompCompMap;
struct malnComp;
struct malnBlk;
struct Genome;
struct Genomes;

/* location of a component in the tree, can also be a bit-set for
 * selection. */
enum mafTreeLoc {
    mafTreeLocUnknown     = 0x0,
    mafTreeLocRoot        = 0x1,
    mafTreeLocInternal    = 0x2,
    mafTreeLocLeaf        = 0x4,
    mafTreeLocAll         = 0xF
};

/* link a tree node with component */
struct mafTreeNodeCompLink {
    int treeOrder;          // DFS post-visit order
    ETree *node;
    struct malnComp *comp;
};


/* Construct a MafTree object from an malnBlk object and a newick tree. */
mafTree *mafTree_constructFromNewick(struct Genomes *genomes, struct malnBlk *blk, char *nhTree);

/* Construct a MafTree object from a malnBlk, assuming the last is the root and all others are descents  */
mafTree *mafTree_constructFromAlignment(struct Genomes *genomes, struct malnBlk *blk, double defaultBranchLength);

/* clone a mafTree */
mafTree *mafTree_clone(mafTree *srcMTree, struct malnBlk *destBlk);

/* destroy a MafTree object */
void mafTree_destruct(mafTree *mTree);

/* join two trees at nodes specified by components */
mafTree *mafTree_join(mafTree *mTree1, struct malnComp *comp1, mafTree *mTree2, struct malnComp *comp2, struct malnCompCompMap *srcDestCompMap);

/* format tree as a newick string */
char *mafTree_format(mafTree *mTree);

/* get string version of a mafTreeLoc */
char *mafTreeLoc_str(enum mafTreeLoc loc);

/* get location type in tree */
enum mafTreeLoc mafTreeNodeCompLink_getLoc(struct mafTreeNodeCompLink *ncLink);

/* Remove a node from the tree and free.  Can't delete the root node. */
void mafTree_deleteNode(mafTree *mTree, struct mafTreeNodeCompLink *ncLink);

/* sort children so tests are reproducible */
void mafTree_sortChildren(mafTree *mTree);

#endif
