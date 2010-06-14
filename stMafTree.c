#include "mafJoinTypes.h"
#include "stMafTree.h"
#include "common.h"
#include "jkmaf.h"
#include "sonLibETree.h"
#include "sonLibList.h"
#include "sonLibString.h"
#include "sonLibHash.h"

/* 
 * treeOrder is the DFS post-traversal order number of the node,
 * which corresponds to the order in the MAF.
 *
 * Each eTree node's client data contains a pointer to the struct nodeCompLink
 * node in the stMafTree.  This is initially copied on clone. and then update.
 */

/*
 * Object for containing and manipulating parsed tree and mapping to the order
 * of the maf rows.  Rows of MAF are in DFS post-traversal order.
 */
struct stMafTree {
    ETree *eTree;
    int numNodes;    // array of nodes, indexed by tree DFS order
    struct nodeCompLink *nodes;
};

/* link a tree node with component coordinates */
struct nodeCompLink {
    int treeOrder;
    ETree *eNode;
    const char *orgSeq;   // not owed here, label from node
    int chromStart;       // component coords
    int chromEnd;
};


/* DFS to fill in table of node links and link back with clientData */
static void fillTreeOrderDFS(stMafTree *mTree, ETree *eNode, int *treeOrder) {
    for (int i = 0; i < eTree_getChildNumber(eNode); i++) {
        fillTreeOrderDFS(mTree, eTree_getChild(eNode, i), treeOrder);
    }
    struct nodeCompLink *ncLink = &(mTree->nodes[*treeOrder]);
    ncLink->treeOrder = *treeOrder;
    ncLink->eNode = eNode;
    ncLink->orgSeq = eTree_getLabel(eNode);
    eTree_setClientData(eNode, ncLink);
    (*treeOrder)++;
}

/* Fill in stMafTree.nodes from tree. */
static void fillTreeOrder(stMafTree *mTree) {
    mTree->numNodes = eTree_getNumNodes(mTree->eTree);
    mTree->nodes = needMem(mTree->numNodes * sizeof(struct nodeCompLink));
    int treeOrder = 0;
    fillTreeOrderDFS(mTree, mTree->eTree, &treeOrder);
}

/* fill in one node link */
static void fillNodeLinkFromMafComp(struct nodeCompLink *ncLink, int treeOrder, struct mafComp *comp) {
    if (!sameString(ncLink->orgSeq, comp->src)) {
        errAbort("MAF node with tree order %d label %s doesn't match component label %s", treeOrder, ncLink->orgSeq, comp->src);
    }
    assert(ncLink->treeOrder == treeOrder);
    ncLink->chromStart = comp->start;
    ncLink->chromEnd = comp->start + comp->size;
    if (comp->strand == '-') {
        reverseIntRange(&ncLink->chromStart, &ncLink->chromEnd, comp->srcSize);
    }
}

/* fill in sequence coordinates from  mafAli */
static void fillNodeLinksFromMafAli(stMafTree *mTree, struct mafAli *ali) {
    if (mTree->numNodes != slCount(ali->components)) {
        errAbort("number of nodes in MAF tree (%d) doesn't match rows in MAF alignment (%d)", mTree->numNodes, slCount(ali));
    }
    struct mafComp *comp = ali->components;
    for (int treeOrder = 0; comp != NULL; treeOrder++, comp = comp->next) {
        fillNodeLinkFromMafComp(&(mTree->nodes[treeOrder]), treeOrder, comp);
    }
}

/* construct a stMafTree from a ETree object */
static stMafTree *stMafTree_construct(ETree *eTree) {
    stMafTree *mTree;
    AllocVar(mTree);
    mTree->eTree = eTree;
    fillTreeOrder(mTree);
    return mTree;
}

/* clone a stMafTree */
stMafTree *stMafTree_clone(stMafTree *srcMTree) {
    stMafTree *mTree = stMafTree_construct(eTree_clone(srcMTree->eTree));
    for (int i = 0; i < srcMTree->numNodes; i++) {
        mTree->nodes[i].chromStart = srcMTree->nodes[i].chromStart;
        mTree->nodes[i].chromEnd = srcMTree->nodes[i].chromEnd;
    }
    return mTree;
}

/* parse a tree from a maf */
static stMafTree *stMafTree_parseFromMaf(struct mafAli *ali) {
    ETree *eTree = eTree_parseNewickString(ali->tree);
    stMafTree *mTree = stMafTree_construct(eTree);
    fillNodeLinksFromMafAli(mTree, ali);
    return mTree;
}

/* construct an ETree node from a mafComp */
static ETree *mafCompToETreeNode(struct mafComp *comp) {
    ETree *eTree = eTree_construct();
    eTree_setLabel(eTree, comp->src);
    return eTree;
}

/* create a tree from a pairwise maf */
static stMafTree *stMafTree_inferFromPairwiseMaf(struct mafAli *ali, double defaultBranchLength) {
    ETree *eRoot = mafCompToETreeNode(ali->components->next);
    ETree *eLeaf = mafCompToETreeNode(ali->components);
    eTree_setParent(eLeaf, eRoot);
    stMafTree *mTree = stMafTree_construct(eRoot);
    fillNodeLinksFromMafAli(mTree, ali);
    eTree_setBranchLength(eLeaf, defaultBranchLength);
    return mTree;
}

stMafTree *stMafTree_constructFromMaf(struct mafAli *ali, double defaultBranchLength) {
    if (ali->tree != NULL) {
        return stMafTree_parseFromMaf(ali);
    } else if (slCount(ali->components) == 2) {
        return stMafTree_inferFromPairwiseMaf(ali, defaultBranchLength);
    } else {
        errAbort("mafAli isn't pairwise and doesn't have a tree");
        return NULL;
    }
}

void stMafTree_destruct(stMafTree *mTree) {
    if (mTree != NULL) {
        eTree_destruct(mTree->eTree);
        free(mTree);
    }
}

/* find a nodeCompLink by component coords */
static struct nodeCompLink *stMafTree_findNodeCompLink(stMafTree *mTree, const char *orgSeq, int chromStart, int chromEnd) {
    for (int i = 0 ; i < mTree->numNodes; i++) {
        struct nodeCompLink *ncLink = &(mTree->nodes[i]);
        if (sameString(ncLink->orgSeq, orgSeq) && (ncLink->chromStart == chromStart) && (ncLink->chromEnd == chromEnd)) {
            return ncLink;
        }
    }
    errAbort("nodeCompLink not found for %s:%d-%d", orgSeq, chromStart, chromEnd);
    return NULL;
}

/* compare function for coordinates by tree order  */
int stMafTree_treeOrderCmp(stMafTree *mTree, const char *orgSeq1, int chromStart1, int chromEnd1, const char *orgSeq2, int chromStart2, int chromEnd2) {
    return stMafTree_findNodeCompLink(mTree, orgSeq1, chromStart1, chromEnd1) - stMafTree_findNodeCompLink(mTree, orgSeq2, chromStart2, chromEnd2);
}

/* validate that a join node is either a root or a leaf */
static void validateJoinNode(ETree *node) {
    if (!((eTree_getParent(node) == NULL) || (eTree_getChildNumber(node) == 0))) {
        errAbort("join node not root or leaf");
    }
}

/* clone a node, recording a mapping of the source nodeCompLink to the new node */
static ETree *cloneNodeRecordLink(ETree *srcENode, stHash *linkMap) {
    ETree *destENode = eTree_cloneNode(srcENode);
    eTree_setClientData(destENode, NULL);
    stHash_insert(linkMap, eTree_getClientData(srcENode), destENode);
    return destENode;
}

/* recursively clone a tree, recording a mapping of the source nodeCompLink to the new node */
static ETree *cloneTreeRecordLink(ETree *node, ETree *parent2, stHash *linkMap) {
    ETree *node2 = cloneNodeRecordLink(node, linkMap);
    eTree_setParent(node2, parent2);
    for (int i = 0; i < eTree_getChildNumber(node); i++) {
        cloneTreeRecordLink(eTree_getChild(node, i), node2, linkMap);
    }
    return node2;
}

/* clone child and append clones to a give parent node */
static void cloneChildren(ETree *srcParent, ETree *destParent, stHash *linkMap) {
    // two linkMap entry will point to join point
    stHash_insert(linkMap, eTree_getClientData(srcParent), destParent);
    for (int i = 0; i < eTree_getChildNumber(srcParent); i++) {
        eTree_setParent(cloneTreeRecordLink(eTree_getChild(srcParent, i), NULL, linkMap), destParent);
    }
}

/* join two trees at a shared root */
static ETree *joinAtRoots(struct nodeCompLink *ncLink1, struct nodeCompLink *ncLink2, stHash *linkMap) {
    ETree *joinedRoot = cloneTreeRecordLink(ncLink1->eNode, NULL, linkMap);
    cloneChildren(ncLink2->eNode, joinedRoot, linkMap);
    return joinedRoot;
}

/* join two trees with a root being attached to a leaf */
static ETree *joinAtLeaf(ETree *root1, struct nodeCompLink *leafNcLink1, struct nodeCompLink *rootNcLink2, stHash *linkMap) {
    ETree *joinedRoot = cloneTreeRecordLink(root1, NULL, linkMap);
    ETree *joinPoint = stHash_search(linkMap, leafNcLink1);
    cloneChildren(rootNcLink2->eNode, joinPoint, linkMap);
    return joinedRoot;
}

/* fill in coordinates in one node link */
static void setJoinNodeLinkCoords(struct nodeCompLink *srcNcLink, stHash *linkMap) {
    ETree *destENode = stHash_search(linkMap, srcNcLink);
    struct nodeCompLink *destNcLink = eTree_getClientData(destENode);
    assert(stString_eq(destNcLink->orgSeq, srcNcLink->orgSeq));
    if (destNcLink->chromStart == destNcLink->chromEnd) {
        destNcLink->chromStart = srcNcLink->chromStart;
        destNcLink->chromEnd = srcNcLink->chromEnd;
    } else {
        destNcLink->chromStart = min(destNcLink->chromStart, srcNcLink->chromStart);
        destNcLink->chromEnd = max(destNcLink->chromEnd, srcNcLink->chromEnd);
    }
}

/* fill in coordinates in node links */
static void fillJoinNodeLinkCoords(stMafTree *srcMTree, stHash *linkMap) {
    for (int i = 0 ; i < srcMTree->numNodes; i++) {
        setJoinNodeLinkCoords(&(srcMTree->nodes[i]), linkMap);
    }
}

/* join two trees at nodes specified by components */
stMafTree *stMafTree_join(stMafTree *mTree1, const char *orgSeq1, int chromStart1, int chromEnd1, stMafTree *mTree2, const char *orgSeq2, int chromStart2, int chromEnd2) {
    if (!stString_eq(orgSeq1, orgSeq2)) {
        errAbort("join nodes are not the same organism.sequence: \"%s\" and \"%s\"", orgSeq1, orgSeq2);
    }
    struct nodeCompLink *ncLink1 = stMafTree_findNodeCompLink(mTree1, orgSeq1, chromStart1, chromEnd1);
    struct nodeCompLink *ncLink2 = stMafTree_findNodeCompLink(mTree2, orgSeq2, chromStart2, chromEnd2);
    validateJoinNode(ncLink1->eNode);
    validateJoinNode(ncLink1->eNode);

    stHash *linkMap = stHash_construct();
    ETree *joinedRoot = NULL;
    if ((eTree_getParent(ncLink1->eNode) == NULL) && (eTree_getParent(ncLink2->eNode) == NULL)) {
        joinedRoot = joinAtRoots(ncLink1, ncLink2, linkMap);
    } else if (eTree_getParent(ncLink1->eNode) == NULL) {
        joinedRoot = joinAtLeaf(mTree2->eTree, ncLink2, ncLink1, linkMap);
    } else if (eTree_getParent(ncLink2->eNode) == NULL) {
        joinedRoot = joinAtLeaf(mTree1->eTree, ncLink1, ncLink2, linkMap);
    } else {
        errAbort("join nodes don't obey rules");
    }
    stMafTree *mTreeJoined = stMafTree_construct(joinedRoot);
    fillJoinNodeLinkCoords(mTree1, linkMap);
    fillJoinNodeLinkCoords(mTree2, linkMap);
    stHash_destruct(linkMap);
    return mTreeJoined;
}

/* format tree as a newick string */
char *stMafTree_format(stMafTree *mTree) {
    return eTree_getNewickTreeString(mTree->eTree);
}
