#include "mafJoinTypes.h"
#include "mafTree.h"
#include "common.h"
#include "genome.h"
#include "malnCompCompMap.h"
#include "malnComp.h"
#include "malnBlk.h"
#include "sonLibETree.h"
#include "sonLibList.h"
#include "sonLibString.h"

/*
 * Object for containing and manipulating parsed tree and mapping to the order
 * of the maf rows.  Rows of MAF are in DFS post-traversal order.
 * Each eTree node's client data contains a pointer to a nodeCompLink
 * object.  This is initially copied on clone and then update.
 * The treeOrder is the DFS post-traversal order number of the node,
 * which corresponds to the order in the MAF.
 */
struct mafTree {
    ETree *tree;
    struct Genomes *genomes;
};

/* constructor, link with node and comp  */
static struct mafTreeNodeCompLink *mafTreeNodeCompLink_construct(int treeOrder, ETree *node, struct malnComp *comp) {
    struct mafTreeNodeCompLink *ncLink;
    AllocVar(ncLink);
    ncLink->treeOrder = treeOrder;
    if (node != NULL) {
        ncLink->node = node;
        eTree_setClientData(node, ncLink);
    }
    if (comp != NULL) {
        ncLink->comp = comp;
        ncLink->comp->ncLink = ncLink;
    }
    return ncLink;
}

/* destructor */
static void mafTreeNodeCompLink_destruct(struct mafTreeNodeCompLink *ncLink) {
    if (ncLink->node != NULL) {
        eTree_setClientData(ncLink->node, NULL);
    }
    if (ncLink->comp != NULL) {
        ncLink->comp->ncLink = NULL;
    }
    freeMem(ncLink);
}

/* DFS to set or check tree order after a join. */
static void setCheckTreeOrderDFS(mafTree *mTree, ETree *node, bool check, int *treeOrder) {
    for (int i = 0; i < eTree_getChildNumber(node); i++) {
        setCheckTreeOrderDFS(mTree, eTree_getChild(node, i), check, treeOrder);
    }
    struct mafTreeNodeCompLink *ncLink = eTree_getClientData(node);
    if (!sameString(ncLink->comp->seq->orgSeqName, eTree_getLabel(node))) {
        errAbort("tree component name \"%s\" doesn't match tree node name \"%s\"", ncLink->comp->seq->orgSeqName, eTree_getLabel(node));
    }
    if (!check) {
        ncLink->treeOrder = *treeOrder;
    } else if (ncLink->treeOrder != *treeOrder) {
        errAbort("expected tree order (%d) doesn't match actual tree node order (%d) for \"%s\"", *treeOrder, ncLink->treeOrder, eTree_getLabel(node));
    }
    (*treeOrder)++;
}
/* Set tree order after in nodeCompLinks after a join. */
static void setCheckTreeOrder(mafTree *mTree, bool check) {
    int treeOrder = 0;
    setCheckTreeOrderDFS(mTree, mTree->tree, check, &treeOrder);
}

/* DFS to fill in table of node links and link back with clientData */
static void fillNodeCompLinksDFS(mafTree *mTree, ETree *node, int *treeOrder, struct malnComp *treeComps[]) {
    for (int i = 0; i < eTree_getChildNumber(node); i++) {
        fillNodeCompLinksDFS(mTree, eTree_getChild(node, i), treeOrder, treeComps);
    }
    struct mafTreeNodeCompLink *ncLink = mafTreeNodeCompLink_construct(*treeOrder, node, treeComps[*treeOrder]);
    (*treeOrder)++;
    if (!sameString(ncLink->comp->seq->orgSeqName, eTree_getLabel(node))) {
        errAbort("tree component name \"%s\" doesn't match tree node name \"%s\"", ncLink->comp->seq->orgSeqName, eTree_getLabel(node));
    }
}

/* Fill in mafTreeNodeCompLinks from tree and components. */
static void fillNodeCompLinks(mafTree *mTree, struct malnBlk *blk) {
    // tree order array of components
    int numNodes = slCount(blk->comps);
    struct malnComp *treeComps[numNodes];
    struct malnComp *comp = blk->comps; 
    // FIXME: don't need array if build in reverse order
    for (int i = 0; i < numNodes; i++, comp = comp->next) {
        treeComps[i] = comp;
    }
    assert(comp == NULL);

    // DFS build of tree
    int treeOrder = 0;
    fillNodeCompLinksDFS(mTree, mTree->tree, &treeOrder, treeComps);
}

/* construct a mafTree from a ETree object */
static mafTree *mafTree_construct(struct Genomes *genomes, ETree *eTree) {
    mafTree *mTree;
    AllocVar(mTree);
    mTree->tree = eTree;
    mTree->genomes = genomes;
    return mTree;
}

/* clone a mafTree */
mafTree *mafTree_clone(mafTree *srcMTree, struct malnBlk *destBlk) {
    mafTree *mTree = mafTree_construct(srcMTree->genomes, eTree_clone(srcMTree->tree));
    fillNodeCompLinks(mTree, destBlk);
    return mTree;
}

/* Construct a MafTree object from an malnBlk object and a newick tree. */
mafTree *mafTree_constructFromNewick(struct Genomes *genomes, struct malnBlk *blk, char *nhTree) {
    mafTree *mTree = mafTree_construct(genomes, eTree_parseNewickString(nhTree));
    fillNodeCompLinks(mTree, blk);
    return mTree;
}

/* construct an ETree node from a malnComp */
static ETree *malnCompToETreeNode(struct malnComp *comp) {
    ETree *eTree = eTree_construct();
    eTree_setLabel(eTree, comp->seq->orgSeqName);
    return eTree;
}

/* Construct a MafTree object from a malnBlk, assuming the last is the root and all others are descents  */
mafTree *mafTree_constructFromAlignment(struct Genomes *genomes, struct malnBlk *blk, double defaultBranchLength) {
    struct malnComp *rootComp = slLastEl(blk->comps);
    ETree *eRoot = malnCompToETreeNode(rootComp);
    for (struct malnComp *comp = blk->comps; comp != rootComp; comp = comp->next) {
        ETree *eLeaf = malnCompToETreeNode(comp);
        eTree_setParent(eLeaf, eRoot);
        eTree_setBranchLength(eLeaf, defaultBranchLength);
    }
    mafTree *mTree = mafTree_construct(genomes, eRoot);
    fillNodeCompLinks(mTree, blk);
    return mTree;
}

/* recursively free  mafTreeNodeCompLink objects */
static void freeMafTreeNodeCompLinks(ETree *node) {
    mafTreeNodeCompLink_destruct(eTree_getClientData(node));
    for (int i = 0; i < eTree_getChildNumber(node); i++) {
        freeMafTreeNodeCompLinks(eTree_getChild(node, i));
    }
}

/* destroy tree and free mafTreeNodeCompLink objects */
void mafTree_destruct(mafTree *mTree) {
    if (mTree != NULL) {
        freeMafTreeNodeCompLinks(mTree->tree);
        eTree_destruct(mTree->tree);
        freeMem(mTree);
    }
}

/* recursively search a tree for node linked to the specified component */
static ETree *searchByComp(ETree *node, struct malnComp *comp) {
    struct mafTreeNodeCompLink *ncLink = eTree_getClientData(node);
    if (ncLink->comp == comp) {
        return node;
    }
    for (int i = 0; i < eTree_getChildNumber(node); i++) {
        ETree *hit = searchByComp(eTree_getChild(node, i), comp);
        if (hit != NULL) {
            return hit;
        }
    }
    return NULL;
}

/* assert nodes are compatible for a join */
static void assertJoinCompatible(ETree *node1, ETree *node2) {
#ifndef NDEBUG
    struct mafTreeNodeCompLink *ncLink1 = eTree_getClientData(node1);
    struct mafTreeNodeCompLink *ncLink2 = eTree_getClientData(node2);
    assert(ncLink1->comp->seq == ncLink2->comp->seq);
#endif
}

/* clone a node, linking in new mafTreeNodeCompLink object */
static ETree *cloneNode(ETree *srcNode, struct malnCompCompMap *srcDestCompMap) {
    ETree *destNode = eTree_cloneNode(srcNode);
    struct mafTreeNodeCompLink *srcNcLink = eTree_getClientData(srcNode);
    struct malnComp *destComp = malnCompCompMap_get(srcDestCompMap, srcNcLink->comp);
    mafTreeNodeCompLink_construct(-1, destNode, destComp);
    return destNode;
}

/* recursively clone a tree */
static ETree *cloneTree(ETree *srcNode, ETree *destParent, struct malnCompCompMap *srcDestCompMap) {
    ETree *destNode = cloneNode(srcNode, srcDestCompMap);
    eTree_setParent(destNode, destParent);
    for (int i = 0; i < eTree_getChildNumber(srcNode); i++) {
        cloneTree(eTree_getChild(srcNode, i), destNode, srcDestCompMap);
    }
    return destNode;
}

/* clone child and append clones to a give parent node */
static void cloneChildren(ETree *srcParent, ETree *destParent, struct malnCompCompMap *srcDestCompMap) {
    for (int i = 0; i < eTree_getChildNumber(srcParent); i++) {
        eTree_setParent(cloneTree(eTree_getChild(srcParent, i), destParent, srcDestCompMap), destParent);
    }
}


/* Clone tree1 and then attach the children of tree2 at the copy of the
 * specified attachment point. */
static ETree *joinTrees(ETree *srcRoot1, ETree *srcAttach1, ETree *srcRoot2, struct malnCompCompMap *srcDestCompMap) {
    assertJoinCompatible(srcAttach1, srcRoot2);
    assert(eTree_getParent(srcRoot2) == NULL);
    ETree *joinedRoot = cloneTree(srcRoot1, NULL, srcDestCompMap);
    struct mafTreeNodeCompLink *srcAttach1NcLink = eTree_getClientData(srcAttach1);
    ETree *joinLeaf = searchByComp(joinedRoot, malnCompCompMap_get(srcDestCompMap, srcAttach1NcLink->comp));
    cloneChildren(srcRoot2, joinLeaf, srcDestCompMap);
    return joinedRoot;
}

/* Join at the specified components, returning new root */
static ETree *joinAtNodes(ETree *root1, ETree *node1, ETree *root2, ETree *node2, struct malnCompCompMap *srcDestCompMap) {
    if ((eTree_getParent(node1) == NULL) && (eTree_getParent(node2) == NULL)) {
        assert((root1 == node1) && (root2 == node2));
        return joinTrees(root1, node1, node2, srcDestCompMap);
    } else if (eTree_getParent(node1) == NULL) {
        assert(root1 == node1);
        return joinTrees(root2, node2, node1, srcDestCompMap);
    } else if (eTree_getParent(node2) == NULL) {
        assert(root2 == node2);
        return joinTrees(root1, node1, node2, srcDestCompMap);
    } else {
        errAbort("join nodes don't obey rules: node1: %s node2: %s", eTree_getLabel(node1), eTree_getLabel(node2));
        return NULL;
    }
}

/* join two trees at nodes specified by components.
 * srcDestCompMap of the src components in both trees to the destination
 * components in the new tree.
 */
mafTree *mafTree_join(mafTree *mTree1, struct malnComp *comp1, mafTree *mTree2, struct malnComp *comp2, struct malnCompCompMap *srcDestCompMap) {
    if (comp1->seq != comp2->seq) {
        errAbort("join nodes are not the same organism.sequence: \"%s\" and \"%s\"", comp1->seq->orgSeqName, comp2->seq->orgSeqName);
    }
    assert(malnCompCompMap_get(srcDestCompMap, comp1) == malnCompCompMap_get(srcDestCompMap, comp2));
    ETree *node1 = searchByComp(mTree1->tree, comp1);
    ETree *node2 = searchByComp(mTree2->tree, comp2);

    ETree *joinedRoot = joinAtNodes(mTree1->tree, node1, mTree2->tree, node2, srcDestCompMap);
    mafTree *mTreeJoined = mafTree_construct(mTree1->genomes, joinedRoot);
    setCheckTreeOrder(mTreeJoined, false);
    return mTreeJoined;
}

/* format tree as a newick string */
char *mafTree_format(mafTree *mTree) {
    return eTree_getNewickTreeString(mTree->tree);
}

/* get string version of a mafTreeLoc */
char *mafTreeLoc_str(enum mafTreeLoc loc) {
    switch (loc) {
    case mafTreeLocUnknown:
        return "unknown";
    case mafTreeLocRoot:
        return "root";
    case mafTreeLocInternal:
        return "internal";
    case mafTreeLocLeaf:
        return "leaf";
    default:
        return "INVALID";
    }
}

/* get location type in tree */
enum mafTreeLoc mafTreeNodeCompLink_getLoc(struct mafTreeNodeCompLink *ncLink) {
    if (eTree_getParent(ncLink->node) == NULL) {
        return mafTreeLocRoot;
    } else if (eTree_getChildNumber(ncLink->node) == 0) {
        return mafTreeLocLeaf;
    } else {
        return mafTreeLocInternal;
    }
}

/* Remove a node from the tree and free.  Can't delete the root node. */
void mafTree_deleteNode(mafTree *mTree, struct mafTreeNodeCompLink *ncLink) {
    ETree *node = ncLink->node;
    ETree *parent = eTree_getParent(node);
    if (parent == NULL) {
        errAbort("BUG: can't remove tree root node");
    }
    eTree_setParent(node, NULL);
    // setParent changes node children
    while (eTree_getChildNumber(node) > 0) {
        eTree_setParent(eTree_getChild(node, 0), parent);
    }
    freeMafTreeNodeCompLinks(node);
    eTree_destruct(node);
    setCheckTreeOrder(mTree, false);
}

/* enforce order of children */
static int sortChildrenCmpFn(ETree *a, ETree *b) {
    int diff = strcmp(eTree_getLabel(a), eTree_getLabel(b));
    if (diff == 0) {
        // same names, sort by seq and location, if available
        struct mafTreeNodeCompLink *ncLinkA = eTree_getClientData(a), *ncLinkB = eTree_getClientData(b);
        if ((ncLinkA != NULL) && (ncLinkB != NULL)) {
            diff = malnComp_cmp(ncLinkA->comp, ncLinkB->comp);
        }
    }
    return diff;
}

/* sort children so tests are reproducible */
void mafTree_sortChildren(mafTree *mTree) {
    eTree_sortChildren(mTree->tree, sortChildrenCmpFn);
    setCheckTreeOrder(mTree, false);
}

/* forward declaration */
static void cloneForSubRangeBlk(ETree *srcParent, ETree *destParent, struct malnCompCompMap *srcDestCompMap);

/* recursive clone a node for subrange block */
static void cloneNodeSubRangeBlk(ETree *srcParent, ETree *srcNode, ETree *destParent, struct malnCompCompMap *srcDestCompMap) {
    struct mafTreeNodeCompLink *srcNcLink = eTree_getClientData(srcNode);
    struct malnComp *destComp = malnCompCompMap_find(srcDestCompMap, srcNcLink->comp);
    if (destComp == NULL) {
        // drop this node
        cloneForSubRangeBlk(srcNode, destParent, srcDestCompMap);
    } else {
        // clone this node
        ETree *destNode = malnCompToETreeNode(destComp);
        eTree_setParent(destNode, destParent);
        eTree_setBranchLength(destNode, eTree_getBranchLength(srcNode));
        mafTreeNodeCompLink_construct(-1, destNode, destComp);
        cloneForSubRangeBlk(srcNode, destNode, srcDestCompMap);
    }
}

/* recursive clone for subrange block */
static void cloneForSubRangeBlk(ETree *srcParent, ETree *destParent, struct malnCompCompMap *srcDestCompMap) {
    for (int i = 0; i < eTree_getChildNumber(srcParent); i++) {
        cloneNodeSubRangeBlk(srcParent, eTree_getChild(srcParent, i), destParent, srcDestCompMap);
    }
}

/*
 * Clone a tree for addition to a subrange block. If if components are in srcDestCompMap
 * clone them, otherwise drop the node and merge the children.
 */
mafTree *mafTree_cloneForSubRangeBlk(mafTree *mTree, struct malnCompCompMap *srcDestCompMap) {
    ETree *destRoot = cloneNode(mTree->tree, srcDestCompMap);
    cloneForSubRangeBlk(mTree->tree, destRoot, srcDestCompMap);
    mafTree *destMTree = mafTree_construct(mTree->genomes,  destRoot);
    setCheckTreeOrder(destMTree, false);
    return destMTree;
}

/* assert sanity of the tree */
void mafTree_assert(mafTree *mTree, struct malnBlk *blk) {
#ifndef NDEBUG
    setCheckTreeOrder(mTree, true);
    assert(eTree_getNumNodes(mTree->tree) == slCount(blk->comps));
#endif
}

/* assert sanity of nodeCompLink */
void mafTreeNodeCompLink_assert(struct mafTreeNodeCompLink *ncLink) {
#ifndef NDEBUG
    if (ncLink != NULL) {
        assert(stString_eq(eTree_getLabel(ncLink->node), ncLink->comp->seq->orgSeqName));
    }
#endif
}
