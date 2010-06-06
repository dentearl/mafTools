#include "stMafTreeTbl.h"
#include "stMafTree.h"
#include "stMafTreeJoined.h"
#include "sonLibList.h"
#include "sonLibETree.h"
#include "common.h"
#include "jkmaf.h"

/* table of trees, used to allocated unique trees, allow for sharing */
struct stMafTreeTbl {
    stList *mTrees;      /* list of stMafTree objects */
    stList *mTreeJoins;  /* list of stMtreeJoined objects */
};

stMafTreeTbl *stMafTreeTbl_construct() {
    stMafTreeTbl *tbl;
    AllocVar(tbl);
    tbl->mTrees = stList_construct();
    return tbl;
}

void stMafTreeTbl_destruct(stMafTreeTbl *tbl) {
    if (tbl != NULL) {
        for (int i = 0; i < stList_length(tbl->mTrees); i++) {
            stMafTree_destruct(stList_get(tbl->mTrees, i));
        }
        for (int i = 0; i < stList_length(tbl->mTreeJoins); i++) {
            stMafTreeJoined_destruct(stList_get(tbl->mTreeJoins, i));
        }
        stList_destruct(tbl->mTrees);
        free(tbl);
    }
}

/* search for a stMafTree, given an eTree */
static stMafTree *findByETree(stMafTreeTbl *tbl, ETree *eTree) {
    for (int i = 0; i < stList_length(tbl->mTrees); i++) {
        stMafTree *mTree = stList_get(tbl->mTrees, i);
        if (eTree_equals(eTree, stMafTree_getETree(mTree))) {
            return mTree;
        }
    }
    return NULL;
}

stMafTree *stMafTreeTbl_obtainForMafAli(stMafTreeTbl *tbl, struct mafAli *ali) {
    if (ali->tree == NULL) {
        errAbort("mafAli does not have a tree");
    }
    ETree *eTree = eTree_parseNewickString(ali->tree);
    stMafTree *mTree = findByETree(tbl, eTree);
    if (mTree == NULL) {
        mTree = stMafTree_constructFromMaf(eTree, ali);
        stList_append(tbl->mTrees, mTree);
    } else {
        eTree_destruct(eTree);
        stMafTree_validateWithMaf(mTree, ali);
    }
    return mTree;
}

/* look for existing joined maftrees */
static stMafTreeJoined *findJoined(stMafTreeTbl *tbl, stMafTree *mTree1, int mafOrder1, stMafTree *mTree2, int mafOrder2) {
    for (int i = 0; i < stList_length(tbl->mTreeJoins); i++) {
        stMafTreeJoined *mtJoined = stList_get(tbl->mTreeJoins, i);
        if (stMafTreeJoined_match(mtJoined, mTree1, mafOrder1, mTree2, mafOrder2)) {
            return mtJoined;
        }
    }
    return NULL;
}

/* add a joined maftrees */
static stMafTreeJoined *addJoined(stMafTreeTbl *tbl, stMafTree *mTree1, int mafOrder1, stMafTree *mTree2, int mafOrder2, stMafTree *mTreeJoined) {
    stMafTreeJoined *mtJoined = stMafTreeJoined_construct(mTree1, mafOrder1, mTree2, mafOrder2, mTreeJoined);
    stList_append(tbl->mTreeJoins, mtJoined);
    return mtJoined;
}

stMafTree *stMafTreeTbl_obtainJoin(stMafTreeTbl *tbl, stMafTree *mTree1, int mafOrder1, stMafTree *mTree2, int mafOrder2) {
    stMafTreeJoined *mtJoined = findJoined(tbl, mTree1, mafOrder1, mTree2, mafOrder2);
    if (mtJoined == NULL) {
        stMafTree *mTreeJoined = stMafTree_join(mTree1, mafOrder1, mTree2, mafOrder2);
        mtJoined = addJoined(tbl, mTree1, mafOrder1, mTree2, mafOrder2, mTreeJoined);
    }
    return stMafTreeJoined_getJoinedMTree(mtJoined);
}
