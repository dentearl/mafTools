#include "mafJoinTypes.h"
#include "stMafTree.h"
#include "common.h"
#include "jkmaf.h"
#include "sonLibETree.h"
#include "sonLibList.h"
#include "sonLibString.h"

/* notes:
 *   - tree designed to be shared by all columns parsed from a given mafAli.
 */

/*
 * Object for containing and manipulating parsed tree and mapping to the order
 * of the maf rows.  Rows of MAF are in DFS post-traversal order.
 */
struct stMafTree {
    ETree *eTree;
};

ETree *stMafTree_getETree(stMafTree *mTree) {
    return mTree->eTree;
}

/* DFS to get MAF order */
static bool findMafOrderDFS(ETree *root, const char *orgSeq, int *curMafOrder) {
    for (int i = 0; i < eTree_getChildNumber(root); i++) {
        if (findMafOrderDFS(eTree_getChild(root, i), orgSeq, curMafOrder)) {
            return true;
        }
    }
    const char *lab = eTree_getLabel(root);
    if ((lab != NULL) && (strcmp(lab, orgSeq) == 0)) {
        return true;
    }
    (*curMafOrder)++;
    return false;
}

/* Get MAF order for an ETree, or -1 if not found. */
static int findMafOrder(ETree *root, const char *orgSeq) {
    int curMafOrder = 0;
    if (findMafOrderDFS(root, orgSeq, &curMafOrder)) {
        return curMafOrder;
    } else {
        return -1;
    }
}

/* Get MAF order for an ETree, or error if not found. */
static int getMafOrder(ETree *root, const char *orgSeq) {
    int mafOrder = findMafOrder(root, orgSeq);
    if (mafOrder < 0) {
        errAbort("%s not found in MAF tree", orgSeq);
    }
    return mafOrder;
}

int stMafTree_findMafOrder(stMafTree *mTree, const char *orgSeq) {
    return findMafOrder(mTree->eTree, orgSeq);
}

int stMafTree_getMafOrder(stMafTree *mTree, const char *orgSeq) {
    return getMafOrder(mTree->eTree, orgSeq);
}

/* Recursively get an eTree node by MAF order number */
static ETree *findByMafOrderDFS(ETree *root, int mafOrder, int *curMafOrder) {
    for (int i = 0; i < eTree_getChildNumber(root); i++) {
        ETree *node = findByMafOrderDFS(eTree_getChild(root, i), mafOrder, curMafOrder);
        if (node != NULL) {
            return node;
        }
    }
    if (*curMafOrder == mafOrder) {
        return root;
    }
    (*curMafOrder)++;
    return NULL;
}

/* Get an eTree node by MAF order number, or NULL if not found */
static ETree *findByMafOrder(ETree *root, int mafOrder) {
    int curMafOrder = 0;
    return findByMafOrderDFS(root, mafOrder, &curMafOrder);
}

/* Recursively get an eTree node by MAF order number, or error if not found */
static ETree *getByMafOrder(ETree *root, int mafOrder) {
    ETree *node = findByMafOrder(root, mafOrder);
    if (node == NULL) {
        errAbort("MAF order %d not found in MAF tree", mafOrder);
    }
    return node;
}

ETree *stMafTree_getByMafOrder(stMafTree *mTree, int mafOrder) {
    return getByMafOrder(mTree->eTree, mafOrder);
}

/* verify a mafComp against MAF order in the tree  */
static void checkCompWithTree(stMafTree *mTree, struct mafComp *comp, int iComp) {
    int mafOrder = stMafTree_getMafOrder(mTree, comp->src);
    if (mafOrder != iComp) {
        errAbort("MAF tree order for %s (%d) doesn't match MAF component order (%d)", comp->src, mafOrder, iComp);
    }
}

/* validate tree MAF order and nodes matches mafAli */
void stMafTree_validateWithMaf(stMafTree *mTree, struct mafAli *ali) {
    if (eTree_getNumNodes(mTree->eTree) != slCount(ali->components)) {
        errAbort("number of nodes in MAF tree (%d) doesn't match rows in MAF alignment (%d)", eTree_getNumNodes(mTree->eTree), slCount(ali));
    }
    struct mafComp *comp = ali->components;
    for (int iComp = 0; comp != NULL; iComp++, comp = comp->next) {
        checkCompWithTree(mTree, comp, iComp);
    }
}

/* construct a stMafTree from a ETree object */
static stMafTree *stMafTree_construct(ETree *eTree) {
    stMafTree *mTree;
    AllocVar(mTree);
    mTree->eTree = eTree;
    return mTree;
}

stMafTree *stMafTree_constructFromMaf(ETree *eTree, struct mafAli *ali) {
    if (ali->tree == NULL) {
        errAbort("mafAli does not have a tree");
    }
    stMafTree *mTree = stMafTree_construct(eTree);
    stMafTree_validateWithMaf(mTree, ali);
    return mTree;
}

void stMafTree_destruct(stMafTree *mTree) {
    if (mTree != NULL) {
        eTree_destruct(mTree->eTree);
        free(mTree);
    }
}

/* validate that a join node is either a root or a leaf */
static void validateJoinNode(ETree *node) {
    if (!((eTree_getParent(node) == NULL) || (eTree_getChildNumber(node) == 0))) {
        errAbort("join node not root or leaf");
    }
}

/* clone child and append clones to a give parent node */
static void cloneChildren(ETree *srcParent, ETree *destParent) {
    for (int i = 0; i < eTree_getChildNumber(srcParent); i++) {
        eTree_setParent(eTree_clone(eTree_getChild(srcParent, i)), destParent);;
    }
}

/* join two trees at a shared root */
static ETree *joinAtRoots(ETree *root1, ETree *root2) {
    ETree *joinedRoot = eTree_clone(root1);
    cloneChildren(root2, joinedRoot);
    return joinedRoot;
}

/* join two trees with a root being attached to a leaf */
static ETree *joinAtLeaf(ETree *root1, int leafMafOrder1, ETree *root2) {
    ETree *joinedRoot = eTree_clone(root1);
    cloneChildren(getByMafOrder(joinedRoot, leafMafOrder1), root2);
    return joinedRoot;
}

/* join two trees at the nodes identified by their maf order */
stMafTree *stMafTree_join(stMafTree *mTree1, int mafOrder1, stMafTree *mTree2, int mafOrder2) {
    ETree *joinNode1 = stMafTree_getByMafOrder(mTree1, mafOrder1);
    ETree *joinNode2 = stMafTree_getByMafOrder(mTree2, mafOrder2);
    validateJoinNode(joinNode1);
    validateJoinNode(joinNode2);
    if (!stString_eq(eTree_getLabel(joinNode1), eTree_getLabel(joinNode2))) {
        errAbort("join nodes are not the same organism.sequence: \"%s\" and \"%s\"", eTree_getLabel(joinNode1), eTree_getLabel(joinNode2));
    }

    ETree *joinedRoot = NULL;
    if ((eTree_getParent(joinNode1) == NULL) && (eTree_getParent(joinNode2) == NULL)) {
        joinedRoot = joinAtRoots(joinNode1, joinNode2);
    } else if (eTree_getParent(joinNode1) == NULL) {
        joinedRoot = joinAtLeaf(mTree2->eTree, mafOrder2, joinNode1);
    } else if (eTree_getParent(joinNode2) == NULL) {
        joinedRoot = joinAtLeaf(mTree1->eTree, mafOrder1, joinNode2);
    } else {
        errAbort("join nodes don't obey rules");
    }
    return stMafTree_construct(joinedRoot);
}
