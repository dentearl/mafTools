#include "stMafTreeJoined.h"
#include "common.h"

/* describes a pair of joined trees.  No objects are owned  */
struct stMafTreeJoined {
    stMafTree *mTree1;
    int mafOrder1;
    stMafTree *mTree2;
    int mafOrder2;
    stMafTree *mTreeJoined;
};

/* constructor  */
stMafTreeJoined *stMafTreeJoined_construct(stMafTree *mTree1, int mafOrder1, stMafTree *mTree2, int mafOrder2, stMafTree *mTreeJoined) {
    stMafTreeJoined *mtJoined;
    AllocVar(mtJoined);
    mtJoined->mTree1 = mTree1;
    mtJoined->mafOrder1 = mafOrder1;
    mtJoined->mTree2 = mTree2;
    mtJoined->mafOrder2 = mafOrder2;
    mtJoined->mTreeJoined = mTreeJoined;
    return mtJoined;
}

/* destructor  */
void stMafTreeJoined_destruct(stMafTreeJoined *mtJoined) {
    if (mtJoined != NULL) {
        free(mtJoined);
    }
}

/* get the joined mTree from an object */
stMafTree *stMafTreeJoined_getJoinedMTree(stMafTreeJoined *mtJoined) {
    return mtJoined->mTreeJoined;
}

/* set the joined mTree in an object */
void stMafTreeJoined_setJoinedMTree(stMafTreeJoined *mtJoined, stMafTree *mTreeJoined) {
    mtJoined->mTreeJoined = mTreeJoined;
}

/* does a  stMafTreeJoined match the key specifications */
bool stMafTreeJoined_match(stMafTreeJoined *mtJoined, stMafTree *mTree1, int mafOrder1, stMafTree *mTree2, int mafOrder2) {
    return (mtJoined->mTree1 == mTree1) && (mtJoined->mafOrder1 == mafOrder1) && (mtJoined->mTree2 == mTree2) && (mtJoined->mafOrder2 == mafOrder2);
}
