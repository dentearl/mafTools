#ifndef stMafTreeJoined_h
#define stMafTreeJoined_h
#include <stdbool.h>
#include "mafJoinTypes.h"

/* constructor, No objects are owned  */
stMafTreeJoined *stMafTreeJoined_construct(stMafTree *mTree1, int mafOrder1, stMafTree *mTree2, int mafOrder2, stMafTree *mTreeJoined);

/* destructor  */
void stMafTreeJoined_destruct(stMafTreeJoined *mtJoined);

/* get the joined mTree from an object */
stMafTree *stMafTreeJoined_getJoinedMTree(stMafTreeJoined *mtJoined);

/* set the joined mTree in an object */
void stMafTreeJoined_setJoinedMTree(stMafTreeJoined *mtJoined, stMafTree *mTreeJoined);

/* does a  stMafTreeJoined match the key specifications */
bool stMafTreeJoined_match(stMafTreeJoined *mtJoined, stMafTree *mTree1, int mafOrder1, stMafTree *mTree2, int mafOrder2);

#endif
